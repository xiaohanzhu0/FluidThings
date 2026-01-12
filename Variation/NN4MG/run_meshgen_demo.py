import argparse
import numpy as np
import torch
import matplotlib.pyplot as plt

from nn4mg import (
    build_boundary_data,
    extract_boundaries,
    init_model_weights,
    load_grid_csv,
    M_fun_torch,
    Model,
    plot_grid,
    train,
)


def main():
    parser = argparse.ArgumentParser(
        description="NN-based variational mesh generation demo."
    )
    parser.add_argument("--x1", default="x1.csv", help="Path to x1.csv")
    parser.add_argument("--x2", default="x2.csv", help="Path to x2.csv")
    parser.add_argument("--steps", type=int, default=40)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--width", type=int, default=128)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--N-bdry", type=int, default=256)
    parser.add_argument("--plot-every", type=int, default=0, help="Plot every N steps during training.")
    parser.add_argument("--no-show", action="store_true", help="Do not show plots.")
    parser.add_argument("--plot-initial", action="store_true", help="Plot initial grid.")
    parser.add_argument("--save-plot", default="", help="Save final grid plot to file.")
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

    train(
        net,
        X1,
        X2,
        boundary,
        M_fun_torch,
        steps=args.steps,
        N_bdry=args.N_bdry,
        lr=args.lr,
        plot_every=args.plot_every or None,
        device=device,
        dtype=dtype,
    )

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
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
