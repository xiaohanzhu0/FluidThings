import numpy as np
import torch

from .boundary import grad_wrt_xy


def compute_misfit_field(
    net,
    X1,
    X2,
    metric_fn=None,
    metric_inv_fn=None,
    Minv_fn=None,
    device=None,
    dtype=None,
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
    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if dtype is None:
        dtype = torch.float32

    xy = np.concatenate(
        (X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1
    )
    xy = torch.from_numpy(xy).to(device=device, dtype=dtype).requires_grad_(True)

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
    else:
        metric_xy = out
        J = None

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
            return (lam_xi**2) * q_xi + (lam_eta**2) * q_eta
        if misfit_type == "target1":
            return (lam_xi**2 * q_xi - 1.0) ** 2 + (lam_eta**2 * q_eta - 1.0) ** 2
        if misfit_type == "target2":
            return torch.abs(lam_xi**2 * q_xi - 1.0) + torch.abs(lam_eta**2 * q_eta - 1.0)
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
    return integrand.detach().cpu().numpy().reshape(X1.shape)


def compute_skewness(X, Y, eps=1e-12):
    dX_dxi = np.zeros_like(X)
    dY_dxi = np.zeros_like(Y)
    dX_deta = np.zeros_like(X)
    dY_deta = np.zeros_like(Y)

    dX_dxi[:, 1:-1] = 0.5 * (X[:, 2:] - X[:, :-2])
    dY_dxi[:, 1:-1] = 0.5 * (Y[:, 2:] - Y[:, :-2])
    dX_dxi[:, 0] = X[:, 1] - X[:, 0]
    dY_dxi[:, 0] = Y[:, 1] - Y[:, 0]
    dX_dxi[:, -1] = X[:, -1] - X[:, -2]
    dY_dxi[:, -1] = Y[:, -1] - Y[:, -2]

    dX_deta[1:-1, :] = 0.5 * (X[2:, :] - X[:-2, :])
    dY_deta[1:-1, :] = 0.5 * (Y[2:, :] - Y[:-2, :])
    dX_deta[0, :] = X[1, :] - X[0, :]
    dY_deta[0, :] = Y[1, :] - Y[0, :]
    dX_deta[-1, :] = X[-1, :] - X[-2, :]
    dY_deta[-1, :] = Y[-1, :] - Y[-2, :]

    v1x, v1y = dX_dxi, dY_dxi
    v2x, v2y = dX_deta, dY_deta
    dot = v1x * v2x + v1y * v2y
    n1 = np.sqrt(v1x * v1x + v1y * v1y)
    n2 = np.sqrt(v2x * v2x + v2y * v2y)
    cosang = dot / (n1 * n2 + eps)
    cosang = np.clip(cosang, -1.0, 1.0)
    angle = np.degrees(np.arccos(cosang))
    skewness = np.abs(90.0 - angle)
    return skewness
