import argparse
from datetime import datetime
from pathlib import Path
import numpy as np
import torch
import matplotlib.pyplot as plt

from nn4mg import (
    compute_misfit_field,
    compute_gradation_field,
    compute_skewness,
    build_boundary_data,
    extract_boundaries,
    init_model_weights,
    load_grid_csv,
    metric_fn_from_npz,
    M_fun_torch,
    M_fun_inv_torch,
    Model,
    plot_grid,
    plot_scalar_field,
    prepare_run_dir,
    save_run,
    train,
)


def main():
    script_dir = Path(__file__).resolve().parent
    default_x1 = script_dir / "x1.csv"
    default_x2 = script_dir / "x2.csv"
    parser = argparse.ArgumentParser(
        description="NN-based variational mesh generation demo."
    )
    parser.add_argument("--x1", default=str(default_x1), help="Path to x1.csv")
    parser.add_argument("--x2", default=str(default_x2), help="Path to x2.csv")
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument(
        "--lr-schedule",
        choices=["constant", "linear", "cosine"],
        default="cosine",
        help="Learning rate schedule.",
    )
    parser.add_argument(
        "--lr-final",
        type=float,
        default=None,
        help="Final learning rate for non-constant schedules (default: lr*0.1).",
    )
    parser.add_argument("--width", type=int, default=64)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--N-int", type=int, default=4096*2)
    parser.add_argument("--N-bdry", type=int, default=256)
    parser.add_argument(
        "--sampling",
        type=int,
        default=1,
        help="Interior sampling mode: 0=grid, 1=uniform unit square.",
    )
    parser.add_argument("--w-neu", type=float, default=0.0)
    parser.add_argument("--w-orth", type=float, default=0.0)
    parser.add_argument("--w-mono", type=float, default=0.0)
    parser.add_argument("--mono-eps", type=float, default=1e-5)
    parser.add_argument("--w-gradation", type=float, default=0.)
    parser.add_argument("--gradation-beta", type=float, default=0.5)
    parser.add_argument("--gradation-beta-weight", type=float, default=1.0)
    parser.add_argument("--det-barrier-scale", type=float, default=1.0)
    parser.add_argument(
        "--misfit",
        choices=["standard", "target1", "target2"],
        default="standard",
        help=(
            "Misfit integrand: "
            "standard=lam^2*q"
            "target1=(lam^2*q-1)^2"
            "target2=abs[(q-mean(q))/mean(q)]"
            "target3=softplus_abs[(q-mean(q))/mean(q)]"
        ),
    )
    parser.add_argument(
        "--problem", type=int, default=2, help="Metric problem id (1 or 2)."
    )
    parser.add_argument(
        "--metric-npz",
        default="./data/harmonic_block9/Mp.npz",
        help="Path to metric npz file (overrides --problem when provided).",
    )
    parser.add_argument(
        "--formulation",
        choices=["primal", "dual"],
        default="primal",
        help="Loss formulation to use (primal: x(s), dual: s(x)).",
    )
    parser.add_argument("--plot-every", type=int, default=0, help="Plot every N steps during training.")
    parser.add_argument("--no-show", action="store_true", help="Do not show plots.")
    parser.add_argument("--plot-initial", action="store_true", help="Plot initial grid.")
    parser.add_argument("--save-plot", default="", help="Save final grid plot to file.")
    parser.add_argument(
        "--save-run-dir",
        default="",
        help="Save model/config to a new subdir under this path.",
    )
    parser.add_argument(
        "--save-run-name",
        default="",
        help="Optional run subdir name (default: timestamp).",
    )
    args = parser.parse_args()

    X1, X2 = load_grid_csv(args.x1, args.x2)
    Nx1, Nx2 = np.shape(X1)

    if args.plot_initial:
        plot_grid(X1, X2, title="Initial grid")
        if not args.no_show:
            plt.show()

    south_np, north_np, west_np, east_np = extract_boundaries(X1, X2)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    dtype = torch.float32
    boundary = build_boundary_data(
        south_np, north_np, west_np, east_np, device=device, dtype=dtype
    )

    net = Model(width=args.width, depth=args.depth)
    init_model_weights(net)

    if args.metric_npz:
        metric_fn_selected = metric_fn_from_npz(args.metric_npz)
        metric_inv_fn_selected = None
    else:
        def metric_fn_selected(x1, x2):
            return M_fun_torch(x1, x2, problem=args.problem)

        def metric_inv_fn_selected(x1, x2):
            return M_fun_inv_torch(x1, x2, problem=args.problem)

    train_params = dict(
        steps=args.steps,
        N_int=args.N_int,
        N_bdry=args.N_bdry,
        lr=args.lr,
        lr_schedule=args.lr_schedule,
        lr_final=args.lr_final,
        plot_every=args.plot_every or None,
        w_neu=args.w_neu,
        w_orth=args.w_orth,
        w_mono=args.w_mono,
        mono_eps=args.mono_eps,
        w_gradation=args.w_gradation,
        gradation_beta=args.gradation_beta,
        gradation_beta_weight=args.gradation_beta_weight,
        det_barrier_scale=args.det_barrier_scale,
        formulation=args.formulation,
        misfit_type=args.misfit,
        sampling=args.sampling,
        device=device,
        dtype=dtype,
    )

    net, best = train(
        net,
        X1,
        X2,
        boundary,
        metric_fn_selected,
        metric_inv_fn=metric_inv_fn_selected,
        **train_params,
    )

    loss_history = best.get("history")
    fig_loss = None
    if loss_history:
        steps = np.asarray(loss_history["step"], dtype=int)
        l_int = np.asarray(loss_history["L_int"], dtype=float)
        l_bdry = np.asarray(loss_history["L_bdry"], dtype=float)
        l_orth = np.asarray(loss_history["L_orth"], dtype=float)
        l_grad = np.asarray(loss_history["L_gradation"], dtype=float)
        l_total = np.asarray(loss_history["L_total"], dtype=float)
        fig_loss, ax_loss = plt.subplots()
        ax_loss.plot(steps, l_total, label="total")
        ax_loss.plot(steps, l_int, label="interior")
        ax_loss.plot(steps, l_bdry, label="boundary")
        ax_loss.plot(steps, l_orth, label="orth")
        ax_loss.plot(steps, l_grad, label="gradation")
        ax_loss.set_xlabel("Step")
        ax_loss.set_ylabel("Loss")
        ax_loss.set_yscale("log")
        ax_loss.legend()
        fig_loss.tight_layout()

    run_dir = None
    if args.save_run_dir:
        run_dir = prepare_run_dir(args.save_run_dir, args.save_run_name or None)
        run_config = {
            "args": vars(args),
            "train": train_params,
            "model": {"width": args.width, "depth": args.depth},
            "best_loss": best["loss"],
            "timestamp": datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        }
        save_run(run_dir, net, run_config)
        if loss_history:
            history_path = run_dir / "loss_history.txt"
            history_table = np.column_stack(
                (steps, l_total, l_int, l_bdry, l_orth, l_grad)
            )
            np.savetxt(
                history_path,
                history_table,
                header="step L_total L_int L_bdry L_orth L_gradation",
                fmt=["%d", "%.8e", "%.8e", "%.8e", "%.8e", "%.8e"],
            )
        if fig_loss is not None:
            fig_loss.savefig(run_dir / "loss_history.png", bbox_inches="tight")
        print(f"Saved run to {run_dir}")

    xy_int = np.concatenate(
        (X1.flatten().reshape(-1, 1), X2.flatten().reshape(-1, 1)), axis=1
    )
    xy_int = torch.from_numpy(xy_int).float().to(device=device, dtype=dtype)

    with torch.no_grad():
        out = net(xy_int)

    X = torch.reshape(out[:, 0], (Nx2, Nx1)).detach().cpu().numpy()
    Y = torch.reshape(out[:, 1], (Nx2, Nx1)).detach().cpu().numpy()

    fig, _ = plot_grid(X, Y, title="Deformed grid")
    if args.save_plot:
        fig.savefig(args.save_plot, bbox_inches="tight")

    misfit_field = compute_misfit_field(
        net,
        X1,
        X2,
        metric_fn_selected,
        metric_inv_fn=metric_inv_fn_selected,
        device=device,
        dtype=dtype,
        formulation=args.formulation,
        misfit_type=args.misfit,
    )
    misfit_avg = float(np.mean(misfit_field))
    misfit_max = float(np.max(misfit_field))
    plot_scalar_field(
        X,
        Y,
        misfit_field,
        title=f"Misfit integrand (avg={misfit_avg:.3e}, max={misfit_max:.3e})",
    )

    skew_field = compute_skewness(X, Y)
    skew_avg = float(np.mean(skew_field))
    skew_max = float(np.max(skew_field))
    plot_scalar_field(
        X,
        Y,
        skew_field,
        title=(
            "Skewness |90-θ| (deg) "
            f"(avg={skew_avg:.2f}, max={skew_max:.2f})"
        ),
        cmap="magma",
    )

    grad_field = compute_gradation_field(
        net,
        X1,
        X2,
        formulation=args.formulation,
        device=device,
        dtype=dtype,
    )
    grad_avg = float(np.mean(grad_field))
    grad_max = float(np.max(grad_field))
    plot_scalar_field(
        X,
        Y,
        grad_field,
        title=(
            "Gradation magnitude |grad log h| "
            f"(avg={grad_avg:.3e}, max={grad_max:.3e})"
        ),
        cmap="inferno",
    )
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
