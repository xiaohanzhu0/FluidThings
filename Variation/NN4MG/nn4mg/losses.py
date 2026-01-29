import torch

from .boundary import grad_wrt_xy


def interior_loss(
    net,
    xy_int,
    metric_fn=None,
    metric_inv_fn=None,
    Minv_fn=None,
    lam_xi=1.0,
    lam_eta=1.0,
    lam_mode="fixed",
    lam_detach=True,
    lam_eps=1e-12,
    formulation="dual",
    eps_det=1e-8,
    w_det_barrier=1.0,
    det_barrier_scale=100.0,
    det_target=1e-6,
    misfit_type="standard",
    w_gradation=0.0,
    gradation_eps=1e-12,
    gradation_beta=0.0,
    gradation_beta_weight=1.0,
):
    if metric_inv_fn is None and Minv_fn is not None:
        metric_inv_fn = Minv_fn
    if metric_fn is None and metric_inv_fn is None:
        raise ValueError("metric_fn or metric_inv_fn must be provided.")
    if formulation not in {"dual", "primal"}:
        raise ValueError(f"Unknown formulation: {formulation}")
    if lam_mode not in {"fixed", "eq9"}:
        raise ValueError(f"Unknown lam_mode: {lam_mode}")

    xy = xy_int.requires_grad_(True)

    out = net(xy)
    xi = out[:, 0:1]
    eta = out[:, 1:2]

    g_xi = grad_wrt_xy(xi, xy)
    g_eta = grad_wrt_xy(eta, xy)

    S11, S12 = g_xi[:, 0], g_xi[:, 1]
    S21, S22 = g_eta[:, 0], g_eta[:, 1]
    detS = S11 * S22 - S12 * S21

    detS_pos = detS.clamp_min(eps_det)
    if formulation == "dual":
        metric_xy = xy
        J = 1.0 / detS_pos
        J_map = J
    else:
        metric_xy = out
        J = 1.0 / detS_pos
        J_map = detS_pos

    x = metric_xy[:, 0]
    y = metric_xy[:, 1]

    def invert_metric_components(M11, M12, M22):
        det = M11 * M22 - M12 * M12
        return M22 / det, -M12 / det, M11 / det

    if metric_fn is not None:
        M11, M12, M22 = metric_fn(x, y)
    else:
        M11_inv, M12_inv, M22_inv = metric_inv_fn(x, y)
        M11, M12, M22 = invert_metric_components(M11_inv, M12_inv, M22_inv)

    if metric_inv_fn is not None:
        M11_inv, M12_inv, M22_inv = metric_inv_fn(x, y)
    else:
        M11_inv, M12_inv, M22_inv = invert_metric_components(M11, M12, M22)

    def quadform_from_components(g, M11, M12, M22):
        gx, gy = g[:, 0], g[:, 1]

        M11 = M11.squeeze(-1)
        M12 = M12.squeeze(-1)
        M22 = M22.squeeze(-1)

        return M11 * gx * gx + 2.0 * M12 * gx * gy + M22 * gy * gy

    def sigma_from_q(q, eps):
        num = q.abs().mean()
        den = (q * q).mean().clamp_min(eps)
        sigma = torch.sqrt(num / den)
        return sigma.detach() if lam_detach else sigma
    
    def misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type):
        if misfit_type == "standard":
            return (lam_xi**2) * q_xi**2 + (lam_eta**2) * q_eta**2
        if misfit_type == "target1":
            return (lam_xi**2 * q_xi - 1.0) ** 2 + (lam_eta**2 * q_eta - 1.0) ** 2
        if misfit_type == "target2":
            return torch.abs((q_xi - q_xi.mean())/q_xi.mean()) + torch.abs((q_eta - q_eta.mean())/q_eta.mean())
        if misfit_type == "target3":
            aux_xi = (q_xi - q_xi.mean())/q_xi.mean()
            aux_eta = (q_eta - q_eta.mean())/q_eta.mean()
            return torch.nn.functional.softplus(aux_xi) + torch.nn.functional.softplus(-aux_xi) + \
                   torch.nn.functional.softplus(aux_eta) + torch.nn.functional.softplus(-aux_eta) - 4*torch.log(torch.tensor(2))
        raise ValueError(f"Unknown misfit_type: {misfit_type}")

    if formulation == "primal":
        col_s1 = torch.cat((g_xi[:, 0:1], g_eta[:, 0:1]), dim=1)
        col_s2 = torch.cat((g_xi[:, 1:2], g_eta[:, 1:2]), dim=1)
        q_xi = quadform_from_components(col_s1, M11, M12, M22)
        q_eta = quadform_from_components(col_s2, M11, M12, M22)
        if lam_mode == "eq9":
            lam_xi = sigma_from_q(q_xi, lam_eps)
            lam_eta = sigma_from_q(q_eta, lam_eps)
        integrand = misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type)
    else:
        q_xi = quadform_from_components(g_xi, M11_inv, M12_inv, M22_inv)
        q_eta = quadform_from_components(g_eta, M11_inv, M12_inv, M22_inv)
        if lam_mode == "eq9":
            lam_xi = sigma_from_q(q_xi, lam_eps)
            lam_eta = sigma_from_q(q_eta, lam_eps)
        integrand = misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type) * J
    loss = integrand.mean()

    L_gradation = torch.zeros((), device=detS.device, dtype=detS.dtype)
    if w_gradation > 0.0:
        log_h = 0.5 * torch.log(J_map.clamp_min(gradation_eps))
        grad_log_h = grad_wrt_xy(log_h, xy)
        grad_mag = torch.sqrt((grad_log_h * grad_log_h).sum(dim=1) + gradation_eps)
        L_gradation = (grad_mag * grad_mag).mean()
        if gradation_beta > 0.0:
            L_gradation = L_gradation + gradation_beta_weight * torch.nn.functional.softplus(
                grad_mag - gradation_beta
            ).pow(2).mean()
        loss = loss + w_gradation * L_gradation

    if w_det_barrier > 0.0:
        loss = loss + det_barrier_scale * w_det_barrier * torch.nn.functional.cross_entropy(
            1 / detS, torch.tensor(1, device=detS.device)
        ).mean()

    _ = J  # preserve original structure from the notebook
    _ = det_target

    stats = {
        "L_int": loss.detach(),
        "L_gradation": L_gradation.detach(),
    }
    return loss, stats


def orth_loss_forward(Fnet, se, eps=1e-12):
    se = se.requires_grad_(True)
    X = Fnet(se)
    x = X[:, 0:1]
    y = X[:, 1:2]

    dx = torch.autograd.grad(x.sum(), se, create_graph=True, retain_graph=True)[0]
    dy = torch.autograd.grad(y.sum(), se, create_graph=True, retain_graph=True)[0]

    F_xi = torch.cat([dx[:, 0:1], dy[:, 0:1]], dim=1)
    F_eta = torch.cat([dx[:, 1:2], dy[:, 1:2]], dim=1)

    g12 = (F_xi * F_eta).sum(dim=1)
    g11 = (F_xi * F_xi).sum(dim=1)
    g22 = (F_eta * F_eta).sum(dim=1)

    return (g12 * g12 / (g11 * g22 + eps)).mean()


def orth_loss_inverse(net, xy, eps=1e-12):
    xy = xy.requires_grad_(True)
    s = net(xy)
    xi = s[:, 0:1]
    eta = s[:, 1:2]

    g_xi = grad_wrt_xy(xi, xy)
    g_eta = grad_wrt_xy(eta, xy)

    g12 = (g_xi * g_eta).sum(dim=1)
    g11 = (g_xi * g_xi).sum(dim=1)
    g22 = (g_eta * g_eta).sum(dim=1)

    return (g12 * g12 / (g11 * g22 + eps)).mean()
