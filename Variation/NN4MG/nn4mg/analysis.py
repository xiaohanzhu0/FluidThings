import numpy as np
import torch

from .boundary import grad_wrt_xy


def compute_misfit_field(
    net, X1, X2, Minv_fn, device=None, dtype=None, lam_xi=1.0, lam_eta=1.0
):
    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if dtype is None:
        dtype = torch.float32

    xy = np.concatenate(
        (X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1
    )
    xy = torch.from_numpy(xy).to(device=device, dtype=dtype).requires_grad_(True)

    s = net(xy)
    xi = s[:, 0:1]
    eta = s[:, 1:2]

    g_xi = grad_wrt_xy(xi, xy)
    g_eta = grad_wrt_xy(eta, xy)

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
