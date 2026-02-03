import math
import numpy as np
import torch

from .losses_decoupled import interior_loss_decoupled


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


def train_decoupled(
    net,
    X1,
    X2,
    metric_fn=None,
    metric_inv_fn=None,
    Minv_fn=None,
    formulation="dual",
    steps=50,
    N_int=4096,
    lr=2e-4,
    lr_schedule="constant",
    lr_final=None,
    lam_xi=1,
    lam_eta=1,
    lam_mode="fixed",
    lam_detach=False,
    lam_eps=1e-12,
    eps_det=1e-8,
    misfit_type="standard",
    grad_clip=None,
    log_every=1,
    device=None,
    dtype=None,
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

    xy_plot = torch.from_numpy(
        np.concatenate((X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1)
    ).float()
    xy_plot = xy_plot.to(device=device, dtype=dtype)

    best = {"loss": float("inf"), "state": None}
    history = {
        "step": [],
        "L_int": [],
        "L_total": [],
    }

    for step in range(1, steps + 1):
        current_lr = lr_at_step(step)
        if schedule_mode != "constant":
            for group in opt.param_groups:
                group["lr"] = current_lr

        if sampling_mode == 0:
            xy_int = xy_plot.detach()
        else:
            xy_int = torch.rand((N_int, 2), device=device, dtype=dtype)

        L_int, _ = interior_loss_decoupled(
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
            eps_det=eps_det,
            misfit_type=misfit_type,
        )

        loss = L_int

        history["step"].append(step)
        history["L_int"].append(L_int.detach().item())
        history["L_total"].append(loss.detach().item())

        opt.zero_grad(set_to_none=True)
        loss.backward()

        if grad_clip is not None:
            torch.nn.utils.clip_grad_norm_(net.parameters(), grad_clip)
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
                f"Lint {L_int.item():.3e}"
            )

    if best["state"] is not None:
        net.load_state_dict(best["state"])
    best["history"] = history
    return net, best
