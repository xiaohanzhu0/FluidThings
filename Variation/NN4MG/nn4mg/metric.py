import numpy as np
import torch


def M_fun(x1, x2):
    M11 = 40000 * (1 + 15 * x1) ** (-2)
    M12 = np.zeros_like(x1)
    M22 = 40000 * (1 + 15 * x2) ** (-2)
    return M11, M12, M22


def M_fun_torch(x1, x2):
    x1_torch = x1
    x2_torch = x2

    M11 = 40000 * (1 + 15 * x1_torch) ** (-2)
    M12 = torch.zeros_like(x1_torch)
    M22 = 40000 * (1 + 15 * x2_torch) ** (-2)
    return M11, M12, M22


def M_fun_inv_torch(x1, x2):
    M11, M12, M22 = M_fun_torch(x1, x2)
    J = M11 * M22 - M12**2
    M11_inv = M22 / J
    M12_inv = -M12 / J
    M22_inv = M11 / J
    return M11_inv, M12_inv, M22_inv
