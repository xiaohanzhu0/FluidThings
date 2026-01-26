#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np


def _normalize_delimiter(raw):
    if raw is None:
        return None
    value = raw.strip()
    if value.lower() in {"auto", "whitespace", "space"}:
        return None
    if value == "\\t":
        return "\t"
    return value


def _load_array(path, delimiter, skiprows, force_2d=True):
    if delimiter is None:
        arr = np.loadtxt(path, skiprows=skiprows)
    else:
        arr = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows)
    if force_2d:
        return np.atleast_2d(arr)
    return arr


def _normalize_vector_or_grid(arr):
    if arr.ndim == 2 and 1 in arr.shape:
        return arr.reshape(-1)
    return arr


def _build_grid_payload(x1_raw, x2_raw):
    x1_raw = _normalize_vector_or_grid(x1_raw)
    x2_raw = _normalize_vector_or_grid(x2_raw)

    if x1_raw.ndim == 1 and x2_raw.ndim == 1:
        X1, X2 = np.meshgrid(x1_raw, x2_raw)
        return {"x1": x1_raw, "x2": x2_raw, "X1": X1, "X2": X2}

    if x1_raw.shape != x2_raw.shape:
        raise ValueError(
            f"Shape mismatch: x1 {x1_raw.shape}, x2 {x2_raw.shape}"
        )
    return {"x1": x1_raw, "x2": x2_raw, "X1": x1_raw, "X2": x2_raw}


def _append_inverse(payload, M11, M12, M22):
    det = M11 * M22 - M12 * M12
    payload.update(
        {
            "det": det,
            "M11_inv": M22 / det,
            "M12_inv": -M12 / det,
            "M22_inv": M11 / det,
        }
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert metric components stored on a uniform unit-square grid "
            "into .npz bundles for NN4MG."
        )
    )
    parser.add_argument("--m11", default="./data/harmonic_transform_output/Mp11.txt", help="Path to Mp11 csv/txt file.")
    parser.add_argument("--m12", default="./data/harmonic_transform_output/Mp12.txt", help="Path to Mp12 csv/txt file.")
    parser.add_argument("--m22", default="./data/harmonic_transform_output/Mp22.txt", help="Path to Mp22 csv/txt file.")
    parser.add_argument("--phys-m11", default="./data/harmonic_transform_output/M11.txt", help="Path to M11 csv/txt file.")
    parser.add_argument("--phys-m12", default="./data/harmonic_transform_output/M12.txt", help="Path to M12 csv/txt file.")
    parser.add_argument("--phys-m22", default="./data/harmonic_transform_output/M22.txt", help="Path to M22 csv/txt file.")
    parser.add_argument("--x1", default="./data/harmonic_transform_output/x1.txt", help="Path to x1 csv/txt file.")
    parser.add_argument("--x2", default="./data/harmonic_transform_output/x2.txt", help="Path to x2 csv/txt file.")
    parser.add_argument(
        "--output",
        default="./data/harmonic_transform_output/Mp.npz",
        help="Output .npz file path for the parametrized-space metric.",
    )
    parser.add_argument(
        "--phys-output",
        default="./data/harmonic_transform_output/M.npz",
        help="Output .npz file path for the physical-space metric.",
    )
    parser.add_argument(
        "--grid-output",
        default="./data/harmonic_transform_output/X.npz",
        help="Output .npz file path for grid locations.",
    )
    parser.add_argument(
        "--delimiter",
        default="auto",
        help="Field delimiter (default: ','). Use 'auto' for whitespace.",
    )
    parser.add_argument(
        "--skiprows",
        type=int,
        default=0,
        help="Rows to skip at the start of each file (default: 0).",
    )
    parser.add_argument(
        "--transpose",
        action="store_true",
        help="Transpose inputs before saving.",
    )
    parser.add_argument(
        "--include-inverse",
        action="store_true",
        help="Include inverse metric components and determinant.",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="Write a compressed .npz file.",
    )
    args = parser.parse_args()

    delimiter = _normalize_delimiter(args.delimiter)
    m11_path = Path(args.m11)
    m12_path = Path(args.m12)
    m22_path = Path(args.m22)
    phys_m11_path = Path(args.phys_m11)
    phys_m12_path = Path(args.phys_m12)
    phys_m22_path = Path(args.phys_m22)
    x1_path = Path(args.x1)
    x2_path = Path(args.x2)

    M11 = _load_array(m11_path, delimiter, args.skiprows, force_2d=True)
    M12 = _load_array(m12_path, delimiter, args.skiprows, force_2d=True)
    M22 = _load_array(m22_path, delimiter, args.skiprows, force_2d=True)
    M11_phys = _load_array(phys_m11_path, delimiter, args.skiprows, force_2d=True)
    M12_phys = _load_array(phys_m12_path, delimiter, args.skiprows, force_2d=True)
    M22_phys = _load_array(phys_m22_path, delimiter, args.skiprows, force_2d=True)
    x1_raw = _load_array(x1_path, delimiter, args.skiprows, force_2d=False)
    x2_raw = _load_array(x2_path, delimiter, args.skiprows, force_2d=False)

    if args.transpose:
        M11 = M11.T
        M12 = M12.T
        M22 = M22.T
        M11_phys = M11_phys.T
        M12_phys = M12_phys.T
        M22_phys = M22_phys.T
        if x1_raw.ndim == 2 and x2_raw.ndim == 2:
            x1_raw = x1_raw.T
            x2_raw = x2_raw.T

    if M11.shape != M12.shape or M11.shape != M22.shape:
        raise ValueError(
            f"Shape mismatch: M11 {M11.shape}, M12 {M12.shape}, M22 {M22.shape}"
        )
    if M11_phys.shape != M12_phys.shape or M11_phys.shape != M22_phys.shape:
        raise ValueError(
            f"Shape mismatch: M11_phys {M11_phys.shape}, "
            f"M12_phys {M12_phys.shape}, M22_phys {M22_phys.shape}"
        )

    n_rows, n_cols = M11.shape
    x = np.linspace(0.0, 1.0, n_cols)
    y = np.linspace(0.0, 1.0, n_rows)
    X1, X2 = np.meshgrid(x, y)

    payload = {
        "x1": x,
        "x2": y,
        "X1": X1,
        "X2": X2,
        "M11": M11,
        "M12": M12,
        "M22": M22,
    }

    if args.include_inverse:
        _append_inverse(payload, M11, M12, M22)

    metric_output = Path(args.output)
    if args.compress:
        np.savez_compressed(metric_output, **payload)
    else:
        np.savez(metric_output, **payload)

    grid_payload = _build_grid_payload(x1_raw, x2_raw)
    grid_shape = grid_payload["X1"].shape
    if M11_phys.shape != grid_shape:
        raise ValueError(
            f"Shape mismatch: M11_phys {M11_phys.shape}, grid {grid_shape}"
        )

    phys_payload = {
        **grid_payload,
        "M11": M11_phys,
        "M12": M12_phys,
        "M22": M22_phys,
    }
    if args.include_inverse:
        _append_inverse(phys_payload, M11_phys, M12_phys, M22_phys)

    phys_output = Path(args.phys_output)
    if args.compress:
        np.savez_compressed(phys_output, **phys_payload)
    else:
        np.savez(phys_output, **phys_payload)

    grid_output = Path(args.grid_output)
    if args.compress:
        np.savez_compressed(grid_output, **grid_payload)
    else:
        np.savez(grid_output, **grid_payload)

    print(
        f"Saved parametrized-space metric to {metric_output} with shape {M11.shape} "
        f"and keys {sorted(payload.keys())}"
    )
    print(
        f"Saved physical-space metric to {phys_output} with shape {M11_phys.shape} "
        f"and keys {sorted(phys_payload.keys())}"
    )
    print(
        f"Saved grid bundle to {grid_output} with keys {sorted(grid_payload.keys())}"
    )


if __name__ == "__main__":
    main()
