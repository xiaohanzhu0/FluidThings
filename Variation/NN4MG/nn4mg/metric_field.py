from typing import Dict, Tuple

import numpy as np
import torch
import torch.nn.functional as F


def _first_available(data: Dict[str, np.ndarray], keys, label):
    for key in keys:
        if key in data:
            return data[key]
    raise KeyError(f"Missing {label} in npz (tried keys: {', '.join(keys)}).")


def _infer_axes(data: Dict[str, np.ndarray], shape: Tuple[int, int]):
    if "x1" in data and "x2" in data:
        x1 = np.asarray(data["x1"]).reshape(-1)
        x2 = np.asarray(data["x2"]).reshape(-1)
        return x1, x2
    if "X1" in data and "X2" in data:
        X1 = np.asarray(data["X1"])
        X2 = np.asarray(data["X2"])
        if X1.shape == shape and X2.shape == shape:
            x1 = X1[0, :]
            x2 = X2[:, 0]
            return np.asarray(x1).reshape(-1), np.asarray(x2).reshape(-1)
    return None, None


def load_metric_npz(path: str) -> Dict[str, np.ndarray]:
    data = np.load(path)
    data_map = {k: data[k] for k in data.files}
    M11 = _first_available(data_map, ("M11", "Mp11"), "M11")
    M12 = _first_available(data_map, ("M12", "Mp12"), "M12")
    M22 = _first_available(data_map, ("M22", "Mp22"), "M22")
    M11 = np.asarray(M11)
    M12 = np.asarray(M12)
    M22 = np.asarray(M22)
    if M11.shape != M12.shape or M11.shape != M22.shape:
        raise ValueError(
            f"Metric component shapes do not match: {M11.shape}, {M12.shape}, {M22.shape}"
        )
    x1, x2 = _infer_axes(data_map, M11.shape)
    return {
        "M11": M11,
        "M12": M12,
        "M22": M22,
        "x1": x1,
        "x2": x2,
    }


def metric_fn_from_npz(path: str):
    payload = load_metric_npz(path)
    M11 = payload["M11"]
    M12 = payload["M12"]
    M22 = payload["M22"]
    x1 = payload["x1"]
    x2 = payload["x2"]

    x_min = float(np.min(x1)) if x1 is not None else 0.0
    x_max = float(np.max(x1)) if x1 is not None else 1.0
    y_min = float(np.min(x2)) if x2 is not None else 0.0
    y_max = float(np.max(x2)) if x2 is not None else 1.0
    if x_max == x_min or y_max == y_min:
        raise ValueError("Invalid axis range for metric field interpolation.")

    grid_np = np.stack([M11, M12, M22], axis=0).astype(np.float32)
    grid_cache = {}

    def metric_fn(x, y):
        if not torch.is_tensor(x) or not torch.is_tensor(y):
            raise TypeError("metric_fn expects torch tensors for x and y.")
        device = x.device
        dtype = x.dtype
        key = (device, dtype)
        if key not in grid_cache:
            grid_cache[key] = torch.from_numpy(grid_np).unsqueeze(0).to(
                device=device, dtype=dtype
            )
        grid = grid_cache[key]

        x_flat = x.reshape(-1)
        y_flat = y.reshape(-1)
        x_norm = (x_flat - x_min) / (x_max - x_min)
        y_norm = (y_flat - y_min) / (y_max - y_min)
        coords = torch.stack((x_norm * 2.0 - 1.0, y_norm * 2.0 - 1.0), dim=-1)
        coords = coords.view(1, -1, 1, 2)
        sample = F.grid_sample(
            grid,
            coords,
            mode="bilinear",
            padding_mode="border",
            align_corners=True,
        )
        sample = sample.squeeze(0).squeeze(-1).transpose(0, 1)
        return sample[:, 0], sample[:, 1], sample[:, 2]

    return metric_fn
