import numpy as np
import torch


def _metric_problem_1_numpy(x1, x2):
    M11 = 40000 * (1 + 15 * x1) ** (-2)
    M12 = np.zeros_like(x1)
    M22 = 40000 * (1 + 15 * x2) ** (-2)
    return M11, M12, M22


def _metric_problem_2_numpy(x1, x2):
    sin_term = np.sin(2 * np.pi * x1) * np.sin(2 * np.pi * x2)
    M11 = 1000 + 600 * sin_term
    M12 = np.zeros_like(x1)
    M22 = 1000 - 600 * sin_term
    return M11, M12, M22


def M_fun(x1, x2, problem=1):
    if problem == 1:
        return _metric_problem_1_numpy(x1, x2)
    if problem == 2:
        return _metric_problem_2_numpy(x1, x2)
    raise ValueError(f"Unknown problem id: {problem}")


def _metric_problem_1_torch(x1, x2):
    x1_torch = x1
    x2_torch = x2

    M11 = 40000 * (1 + 15 * x1_torch) ** (-2)
    M12 = torch.zeros_like(x1_torch)
    M22 = 40000 * (1 + 15 * x2_torch) ** (-2)
    return M11, M12, M22


def _metric_problem_2_torch(x1, x2):
    x1_torch = x1
    x2_torch = x2

    sin_term = torch.sin(2 * torch.pi * x1_torch) * torch.sin(
        2 * torch.pi * x2_torch
    )
    M11 = 1000 + 600 * sin_term
    M12 = torch.zeros_like(x1_torch)
    M22 = 1000 - 600 * sin_term
    return M11, M12, M22


def M_fun_torch(x1, x2, problem=1):
    if problem == 1:
        return _metric_problem_1_torch(x1, x2)
    if problem == 2:
        return _metric_problem_2_torch(x1, x2)
    raise ValueError(f"Unknown problem id: {problem}")


def M_fun_inv_torch(x1, x2, problem=1):
    M11, M12, M22 = M_fun_torch(x1, x2, problem=problem)
    J = M11 * M22 - M12**2
    M11_inv = M22 / J
    M12_inv = -M12 / J
    M22_inv = M11 / J
    return M11_inv, M12_inv, M22_inv
