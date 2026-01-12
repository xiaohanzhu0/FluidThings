import torch

from .boundary import grad_wrt_xy


def interior_loss(
    net,
    xy_int,
    Minv_fn,
    lam_xi=1.0,
    lam_eta=1.0,
    eps_det=1e-8,
    w_det_barrier=1.0,
    det_target=1e-6,
):
    xy = xy_int.requires_grad_(True)

    s = net(xy)

    xi = s[:, 0:1]
    eta = s[:, 1:2]

    g_xi = grad_wrt_xy(xi, xy)
    g_eta = grad_wrt_xy(eta, xy)

    S11, S12 = g_xi[:, 0], g_xi[:, 1]
    S21, S22 = g_eta[:, 0], g_eta[:, 1]
    detS = S11 * S22 - S12 * S21

    detS_pos = detS.clamp_min(eps_det)
    J = 1.0 / detS_pos

    x = xy[:, 0]
    y = xy[:, 1]
    M11_inv, M12_inv, M22_inv = Minv_fn(x, y)

    def quadform_from_components(g, M11_inv, M12_inv, M22_inv):
        gx, gy = g[:, 0], g[:, 1]

        M11 = M11_inv.squeeze(-1)
        M12 = M12_inv.squeeze(-1)
        M22 = M22_inv.squeeze(-1)

        return M11 * gx * gx + 2.0 * M12 * gx * gy + M22 * gy * gy

    q_xi = quadform_from_components(
        torch.cat((g_xi[:, 0:1], g_eta[:, 0:1]), dim=1),
        M11_inv,
        M12_inv,
        M22_inv,
    )
    q_eta = quadform_from_components(
        torch.cat((g_xi[:, 1:2], g_eta[:, 1:2]), dim=1),
        M11_inv,
        M12_inv,
        M22_inv,
    )

    integrand = (lam_xi**2) * q_xi + (lam_eta**2) * q_eta
    loss = integrand.mean()

    if w_det_barrier > 0.0:
        loss = loss + 1 * w_det_barrier * torch.nn.functional.cross_entropy(
            1 / detS, torch.tensor(1, device=detS.device)
        ).mean()

    _ = J  # preserve original structure from the notebook
    _ = det_target

    return loss, {"L_int": loss.detach()}


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
