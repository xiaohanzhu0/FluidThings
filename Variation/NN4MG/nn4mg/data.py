import numpy as np


def load_grid_csv(x1_path: str, x2_path: str):
    X1 = np.loadtxt(x1_path, delimiter=",")
    X2 = np.loadtxt(x2_path, delimiter=",")
    if X1.shape != X2.shape:
        raise ValueError(f"X1/X2 shape mismatch: {X1.shape} vs {X2.shape}")
    return X1, X2


def extract_boundaries(X1: np.ndarray, X2: np.ndarray):
    Nx1, Nx2 = np.shape(X1)
    south_np = np.concatenate((X1[0:1, :], X2[0:1, :]), axis=0).T
    north_np = np.concatenate((X1[-1:Nx2, :], X2[-1:Nx2, :]), axis=0).T
    west_np = np.concatenate((X1[:, 0:1], X2[:, 0:1]), axis=1)
    east_np = np.concatenate((X1[:, -1:Nx1], X2[:, -1:Nx1]), axis=1)
    return south_np, north_np, west_np, east_np
