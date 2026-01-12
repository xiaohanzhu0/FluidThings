import numpy as np
import torch

from .boundary import boundary_loss
from .losses import interior_loss, orth_loss_forward


def train(
    net,
    X1,
    X2,
    boundary,
    Minv_fn,
    steps=50,
    N_int=4096,
    N_bdry=256,
    lr=2e-4,
    lam_xi=1.0,
    lam_eta=1.0,
    w_dir=10.0,
    w_neu=1.0,
    w_det_init=0.1,
    w_det_final=10.0,
    det_target_init=1e-10,
    det_target_final=1e-6,
    grad_clip=1.0,
    log_every=1,
    plot_every=None,
    device=None,
    dtype=None,
):
    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if dtype is None:
        dtype = torch.float32

    net.to(device=device, dtype=dtype)
    opt = torch.optim.Adam(net.parameters(), lr=lr)

    best = {"loss": float("inf"), "state": None}

    fig = None
    ax = None
    for step in range(1, steps + 1):
        xy_int = torch.from_numpy(
            np.concatenate((X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1)
        ).float()
        xy_int = xy_int.to(device=device, dtype=dtype)

        L_int, stats_int = interior_loss(
            net,
            xy_int,
            Minv_fn,
            lam_xi=lam_xi,
            lam_eta=lam_eta,
        )

        L_bdry, stats_bdry = boundary_loss(
            net, boundary, Nb=N_bdry, w_dir=w_dir, w_neu=w_neu
        )

        loss = L_int + L_bdry

        xyS = np.concatenate(
            (X1[0:1, :].flatten().reshape(-1, 1), X2[0:1, :].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyN = np.concatenate(
            (X1[99:100, :].flatten().reshape(-1, 1), X2[99:100, :].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyW = np.concatenate(
            (X1[:, 0:1].flatten().reshape(-1, 1), X2[:, 0:1].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyE = np.concatenate(
            (X1[:, 99:100].flatten().reshape(-1, 1), X2[:, 99:100].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyS = torch.from_numpy(xyS).float().to(device=device, dtype=dtype).requires_grad_(True)
        xyN = torch.from_numpy(xyN).float().to(device=device, dtype=dtype).requires_grad_(True)
        xyW = torch.from_numpy(xyW).float().to(device=device, dtype=dtype).requires_grad_(True)
        xyE = torch.from_numpy(xyE).float().to(device=device, dtype=dtype).requires_grad_(True)

        xiS = net(xyS)[:, 0:1]
        xiN = net(xyN)[:, 0:1]
        etaW = net(xyW)[:, 1:2]
        etaE = net(xyE)[:, 1:2]
        mono_eps = torch.tensor(0.00001, device=device, dtype=dtype)

        d_xi_dx_S = torch.autograd.grad(xiS.sum(), xyS, create_graph=True)[0][:, 0:1]
        d_xi_dx_N = torch.autograd.grad(xiN.sum(), xyN, create_graph=True)[0][:, 0:1]
        d_eta_dy_W = torch.autograd.grad(etaW.sum(), xyW, create_graph=True)[0][:, 1:2]
        d_eta_dy_E = torch.autograd.grad(etaE.sum(), xyE, create_graph=True)[0][:, 1:2]

        L_mono = (
            torch.nn.functional.softplus(mono_eps - d_xi_dx_S).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps - d_xi_dx_N).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps - d_eta_dy_W).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps - d_eta_dy_E).pow(2).mean()
        )
        loss = loss + 0 * L_mono

        w_orth = 1.0
        L_orth = orth_loss_forward(net, xy_int)
        loss = loss + w_orth * L_orth

        opt.zero_grad(set_to_none=True)
        loss.backward()

        if grad_clip is not None:
            pass

        opt.step()

        if loss.item() < best["loss"]:
            best["loss"] = loss.item()
            best["state"] = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}

        if step % log_every == 0:
            print(
                f"step {step:6d} | "
                f"loss {loss.item():.4e} | "
                f"Lint {stats_int['L_int'].item():.3e} "
                f"Lbdry {(L_bdry.detach()).item():.3e} | "
            )
        if plot_every and step % plot_every == 0:
            from .plot import plot_grid
            import matplotlib.pyplot as plt

            Nx1, Nx2 = np.shape(X1)
            with torch.no_grad():
                out = net(xy_int)
            X = torch.reshape(out[:, 0], (Nx2, Nx1)).detach().cpu().numpy()
            Y = torch.reshape(out[:, 1], (Nx2, Nx1)).detach().cpu().numpy()
            if fig is None or ax is None:
                plt.ioff()
                fig, ax = plot_grid(X, Y, title=f"Step {step}")
            else:
                ax.clear()
                ax.plot(X, Y, "k", linewidth=0.1)
                ax.plot(X.T, Y.T, "k", linewidth=0.1)
                ax.set_aspect("equal", adjustable="box")
                ax.set_title(f"Step {step}")

    if best["state"] is not None:
        net.load_state_dict(best["state"])
    return net, best
