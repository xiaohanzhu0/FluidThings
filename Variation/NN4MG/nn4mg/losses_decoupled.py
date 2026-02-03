import math
import torch
import time


def _invert_metric_components(M11, M12, M22):
    det = M11 * M22 - M12 * M12
    return M22 / det, -M12 / det, M11 / det


def _sigma_from_q(q, eps, detach):
    num = q.abs().mean()
    den = (q * q).mean().clamp_min(eps)
    sigma = torch.sqrt(num / den)
    return sigma.detach() if detach else sigma


def _misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type):
    if misfit_type == "standard":
        return (lam_xi**2) * q_xi**2 + (lam_eta**2) * q_eta**2
    if misfit_type == "target1":
        return (lam_xi**2 * q_xi - 1.0) ** 2 + (lam_eta**2 * q_eta - 1.0) ** 2
    if misfit_type == "target2":
        return torch.abs((q_xi - q_xi.mean()) / q_xi.mean()) + torch.abs(
            (q_eta - q_eta.mean()) / q_eta.mean()
        )
    if misfit_type == "target3":
        aux_xi = (q_xi - q_xi.mean()) / q_xi.mean()
        aux_eta = (q_eta - q_eta.mean()) / q_eta.mean()
        return (
            torch.nn.functional.softplus(aux_xi)
            + torch.nn.functional.softplus(-aux_xi)
            + torch.nn.functional.softplus(aux_eta)
            + torch.nn.functional.softplus(-aux_eta)
            - 4.0 * math.log(2.0)
        )
    raise ValueError(f"Unknown misfit_type: {misfit_type}")


def interior_loss_decoupled(
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
    misfit_type="standard",
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

    t0 = time.perf_counter()
    out = net(xy)

    xi = out[:, 0:1]
    eta = out[:, 1:2]

    dxi = torch.autograd.grad(xi.sum(), xy, create_graph=True, retain_graph=True)[0]
    deta = torch.autograd.grad(eta.sum(), xy, create_graph=True, retain_graph=True)[0]
    dxi_dx = dxi[:, 0]
    deta_dy = deta[:, 1]

    dt = time.perf_counter() - t0
    print(f"forward: {dt:.6f} s")

    if formulation == "dual":
        metric_xy = xy
        if metric_inv_fn is not None:
            M11_inv, _, M22_inv = metric_inv_fn(metric_xy[:, 0], metric_xy[:, 1])
        else:
            M11, M12, M22 = metric_fn(metric_xy[:, 0], metric_xy[:, 1])
            M11_inv, _, M22_inv = _invert_metric_components(M11, M12, M22)

        q_xi = M11_inv.squeeze(-1) * dxi_dx * dxi_dx
        q_eta = M22_inv.squeeze(-1) * deta_dy * deta_dy

        if lam_mode == "eq9":
            lam_xi = _sigma_from_q(q_xi, lam_eps, lam_detach)
            lam_eta = _sigma_from_q(q_eta, lam_eps, lam_detach)

        detS = dxi_dx * deta_dy
        J = 1.0 / detS.clamp_min(eps_det)
        integrand = _misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type) * J
    else:
        metric_xy = out
        if metric_fn is not None:
            M11, M12, M22 = metric_fn(metric_xy[:, 0], metric_xy[:, 1])
        else:
            M11_inv, M12_inv, M22_inv = metric_inv_fn(metric_xy[:, 0], metric_xy[:, 1])
            M11, _, M22 = _invert_metric_components(M11_inv, M12_inv, M22_inv)

        q_xi = M11.squeeze(-1) * dxi_dx * dxi_dx
        q_eta = M22.squeeze(-1) * deta_dy * deta_dy

        if lam_mode == "eq9":
            lam_xi = _sigma_from_q(q_xi, lam_eps, lam_detach)
            lam_eta = _sigma_from_q(q_eta, lam_eps, lam_detach)

        integrand = _misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type)

    loss = integrand.mean()
    stats = {"L_int": loss.detach()}
    return loss, stats
