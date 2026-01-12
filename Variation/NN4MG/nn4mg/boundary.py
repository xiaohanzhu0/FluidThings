from dataclasses import dataclass
import torch


@dataclass
class BoundaryData:
    south_xy: torch.Tensor
    north_xy: torch.Tensor
    west_xy: torch.Tensor
    east_xy: torch.Tensor
    n_south: torch.Tensor
    n_north: torch.Tensor
    n_west: torch.Tensor
    n_east: torch.Tensor


def to_xy_tensor(arr, device, dtype, N=None):
    a = arr
    if getattr(a, "ndim", None) is None:
        a = a.numpy()
    if a.ndim == 1:
        if N is None:
            raise AssertionError("Need N to reshape a flat (2*N,) boundary array")
        a = a.reshape(N, 2)
    if a.shape[1] != 2:
        raise AssertionError("Boundary array must have shape (N,2)")
    return torch.as_tensor(a, device=device, dtype=dtype)


def polyline_normals(xy, interior_point):
    d = torch.zeros_like(xy)
    d[1:-1] = xy[2:] - xy[:-2]
    d[0] = xy[1] - xy[0]
    d[-1] = xy[-1] - xy[-2]

    t = d / (d.norm(dim=1, keepdim=True) + 1e-12)
    n = torch.stack([t[:, 1], -t[:, 0]], dim=1)

    v = (interior_point[None, :] - xy)
    flip = (n * v).sum(dim=1, keepdim=True) > 0
    n = torch.where(flip, -n, n)

    n_hat = n / (n.norm(dim=1, keepdim=True) + 1e-12)
    return n_hat


def build_boundary_data(south_np, north_np, west_np, east_np, device, dtype):
    south_xy = to_xy_tensor(south_np, device=device, dtype=dtype)
    north_xy = to_xy_tensor(north_np, device=device, dtype=dtype)
    west_xy = to_xy_tensor(west_np, device=device, dtype=dtype)
    east_xy = to_xy_tensor(east_np, device=device, dtype=dtype)

    all_bdry = torch.cat([south_xy, north_xy, west_xy, east_xy], dim=0)
    interior = all_bdry.mean(dim=0)
    n_south = polyline_normals(south_xy, interior)
    n_north = polyline_normals(north_xy, interior)
    n_west = polyline_normals(west_xy, interior)
    n_east = polyline_normals(east_xy, interior)

    return BoundaryData(
        south_xy=south_xy,
        north_xy=north_xy,
        west_xy=west_xy,
        east_xy=east_xy,
        n_south=n_south,
        n_north=n_north,
        n_west=n_west,
        n_east=n_east,
    )


def sample_boundary(xy, n_hat, Nb):
    idx = torch.randint(0, xy.shape[0], (Nb,), device=xy.device)
    return xy[idx], n_hat[idx]


def grad_wrt_xy(scalar, xy):
    g = torch.autograd.grad(
        outputs=scalar.sum(),
        inputs=xy,
        create_graph=True,
        retain_graph=True,
    )[0]
    return g


def normal_derivative(scalar, xy, n_hat):
    g = grad_wrt_xy(scalar, xy)
    return (g * n_hat).sum(dim=1, keepdim=True)


def boundary_loss(net, boundary: BoundaryData, Nb=256, w_dir=10.0, w_neu=1.0):
    xyS, nS = sample_boundary(boundary.south_xy, boundary.n_south, Nb)
    xyN, nN = sample_boundary(boundary.north_xy, boundary.n_north, Nb)
    xyW, nW = sample_boundary(boundary.west_xy, boundary.n_west, Nb)
    xyE, nE = sample_boundary(boundary.east_xy, boundary.n_east, Nb)

    xyS = xyS.requires_grad_(True)
    xyN = xyN.requires_grad_(True)
    xyW = xyW.requires_grad_(True)
    xyE = xyE.requires_grad_(True)

    xiS, etaS = net(xyS).split(1, dim=1)
    xiN, etaN = net(xyN).split(1, dim=1)
    xiW, etaW = net(xyW).split(1, dim=1)
    xiE, etaE = net(xyE).split(1, dim=1)

    L_dir = (
        torch.nn.functional.mse_loss(etaS, torch.zeros_like(etaS))
        + torch.nn.functional.mse_loss(etaN, torch.ones_like(etaN))
        + torch.nn.functional.mse_loss(xiW, torch.zeros_like(xiW))
        + torch.nn.functional.mse_loss(xiE, torch.ones_like(xiE))
    )

    dxi_dn_S = normal_derivative(xiS, xyS, nS)
    dxi_dn_N = normal_derivative(xiN, xyN, nN)
    deta_dn_W = normal_derivative(etaW, xyW, nW)
    deta_dn_E = normal_derivative(etaE, xyE, nE)

    L_neu = (
        (dxi_dn_S**2).mean()
        + (dxi_dn_N**2).mean()
        + (deta_dn_W**2).mean()
        + (deta_dn_E**2).mean()
    )

    return 0 * w_dir * L_dir + 1000000 * w_neu * L_neu, {
        "L_dir": L_dir.detach(),
        "L_neu": L_neu.detach(),
    }
