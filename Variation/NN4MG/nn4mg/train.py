import math
import numpy as np
import torch

from .boundary import boundary_loss
from .losses import interior_loss, orth_loss_forward, orth_loss_inverse


def _resolve_sampling_mode(sampling):
    if not isinstance(sampling, int):
        raise TypeError("sampling must be an int.")
    if sampling not in (0, 1):
        raise ValueError("Unknown sampling option. Use 0 (grid) or 1 (uniform).")
    return sampling


def _resolve_lr_schedule(lr_schedule):
    if lr_schedule is None:
        return "constant"
    if not isinstance(lr_schedule, str):
        raise TypeError("lr_schedule must be a string.")
    lr_schedule = lr_schedule.lower()
    if lr_schedule not in {"constant", "linear", "cosine"}:
        raise ValueError("Unknown lr_schedule. Use constant, linear, or cosine.")
    return lr_schedule


def train(
    net,
    X1,
    X2,
    boundary,
    metric_fn=None,
    metric_inv_fn=None,
    Minv_fn=None,
    formulation="dual",
    steps=50,
    N_int=4096,
    N_bdry=256,
    lr=2e-4,
    lr_schedule="constant",
    lr_final=None,
    lam_xi=1,
    lam_eta=1,
    lam_mode="fixed",
    lam_detach=False,
    lam_eps=1e-12,
    w_neu=100.0,
    w_det_init=0.1,
    w_det_final=10.0,
    det_target_init=1e-10,
    det_target_final=1e-6,
    w_orth=1.0,
    w_mono=0.0,
    mono_eps=1e-5,
    w_gradation=0.0,
    gradation_eps=1e-12,
    gradation_beta=0.0,
    gradation_beta_weight=1.0,
    det_barrier_scale=100.0,
    grad_clip=1.0,
    log_every=10,
    plot_every=None,
    device=None,
    dtype=None,
    misfit_type="standard",
    sampling=0,
):
    if metric_inv_fn is None and Minv_fn is not None:
        metric_inv_fn = Minv_fn
    if metric_fn is None and metric_inv_fn is None:
        raise ValueError("metric_fn or metric_inv_fn must be provided.")
    if formulation not in {"dual", "primal"}:
        raise ValueError(f"Unknown formulation: {formulation}")

    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if dtype is None:
        dtype = torch.float32

    sampling_mode = _resolve_sampling_mode(sampling)
    if sampling_mode == 1 and (N_int is None or N_int <= 0):
        raise ValueError("N_int must be positive for uniform sampling.")

    net.to(device=device, dtype=dtype)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    schedule_mode = _resolve_lr_schedule(lr_schedule)
    if schedule_mode != "constant":
        if lr_final is None:
            lr_final = lr * 0.1
        lr_final = float(lr_final)
        if lr_final < 0:
            raise ValueError("lr_final must be non-negative.")

    def lr_at_step(step):
        if schedule_mode == "constant":
            return lr
        if steps <= 1:
            return lr_final
        progress = (step - 1) / (steps - 1)
        if schedule_mode == "linear":
            return lr + (lr_final - lr) * progress
        return lr_final + (lr - lr_final) * (0.5 * (1.0 + math.cos(math.pi * progress)))

    best = {"loss": float("inf"), "state": None}
    history = {
        "step": [],
        "L_int": [],
        "L_bdry": [],
        "L_orth": [],
        "L_gradation": [],
        "L_total": [],
    }

    xy_plot = torch.from_numpy(
        np.concatenate((X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1)
    ).float()
    xy_plot = xy_plot.to(device=device, dtype=dtype)

    fig = None
    ax = None
    for step in range(1, steps + 1):
        current_lr = lr_at_step(step)
        if schedule_mode != "constant":
            for group in opt.param_groups:
                group["lr"] = current_lr

        if sampling_mode == 0:
            xy_int = xy_plot.detach()
        else:
            xy_int = torch.rand((N_int, 2), device=device, dtype=dtype)

        L_int, stats_int = interior_loss(
            net,
            xy_int,
            metric_fn=metric_fn,
            metric_inv_fn=metric_inv_fn,
            lam_xi=lam_xi,
            lam_eta=lam_eta,
            lam_mode=lam_mode,
            lam_detach=lam_detach,
            lam_eps=lam_eps,
            formulation=formulation,
            det_barrier_scale=det_barrier_scale,
            misfit_type=misfit_type,
            w_gradation=w_gradation,
            gradation_eps=gradation_eps,
            gradation_beta=gradation_beta,
            gradation_beta_weight=gradation_beta_weight,
        )

        L_bdry, stats_bdry = boundary_loss(net, boundary, Nb=N_bdry, w_neu=w_neu)

        loss = L_int + L_bdry

        n_rows, n_cols = X1.shape
        last_row = n_rows - 1
        last_col = n_cols - 1
        xyS = np.concatenate(
            (X1[0:1, :].flatten().reshape(-1, 1), X2[0:1, :].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyN = np.concatenate(
            (
                X1[last_row : last_row + 1, :].flatten().reshape(-1, 1),
                X2[last_row : last_row + 1, :].flatten().reshape(-1, 1),
            ),
            axis=1,
        )
        xyW = np.concatenate(
            (X1[:, 0:1].flatten().reshape(-1, 1), X2[:, 0:1].flatten().reshape(-1, 1)),
            axis=1,
        )
        xyE = np.concatenate(
            (
                X1[:, last_col : last_col + 1].flatten().reshape(-1, 1),
                X2[:, last_col : last_col + 1].flatten().reshape(-1, 1),
            ),
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
        mono_eps_t = torch.tensor(mono_eps, device=device, dtype=dtype)

        d_xi_dx_S = torch.autograd.grad(xiS.sum(), xyS, create_graph=True)[0][:, 0:1]
        d_xi_dx_N = torch.autograd.grad(xiN.sum(), xyN, create_graph=True)[0][:, 0:1]
        d_eta_dy_W = torch.autograd.grad(etaW.sum(), xyW, create_graph=True)[0][:, 1:2]
        d_eta_dy_E = torch.autograd.grad(etaE.sum(), xyE, create_graph=True)[0][:, 1:2]

        L_mono = (
            torch.nn.functional.softplus(mono_eps_t - d_xi_dx_S).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps_t - d_xi_dx_N).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps_t - d_eta_dy_W).pow(2).mean()
            + torch.nn.functional.softplus(mono_eps_t - d_eta_dy_E).pow(2).mean()
        )
        loss = loss + w_mono * L_mono

        if formulation == "dual":
            L_orth = orth_loss_inverse(net, xy_int)
        else:
            L_orth = orth_loss_forward(net, xy_int)
        loss = loss + w_orth * L_orth
        L_int_val = stats_int["L_int"].item()
        L_grad_val = (w_gradation * stats_int["L_gradation"]).item()
        L_bdry_val = L_bdry.detach().item()
        L_orth_val = (w_orth * L_orth).detach().item()
        history["step"].append(step)
        history["L_int"].append(L_int_val)
        history["L_bdry"].append(L_bdry_val)
        history["L_orth"].append(L_orth_val)
        history["L_gradation"].append(L_grad_val)
        history["L_total"].append(loss.item())

        opt.zero_grad(set_to_none=True)
        loss.backward()

        if grad_clip is not None:
            pass

        opt.step()

        if loss.item() < best["loss"]:
            best["loss"] = loss.item()
            best["state"] = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}

        if step % log_every == 0:
            lr_msg = f"lr {current_lr:.2e} | " if schedule_mode != "constant" else ""
            print(
                f"step {step:6d} | "
                f"{lr_msg}"
                f"loss {loss.item():.4e} | "
                f"Lint {L_int_val:.3e} "
                f"Lbdry {L_bdry_val:.3e} "
                f"Lorth {L_orth_val:.3e} "
                f"Lgrad {L_grad_val:.3e} | "
            )
        if plot_every and step % plot_every == 0:
            from .plot import plot_grid
            import matplotlib.pyplot as plt

            Nx1, Nx2 = np.shape(X1)
            with torch.no_grad():
                out = net(xy_plot)
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
    best["history"] = history
    return net, best
