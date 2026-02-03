import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class Monotonic1D(nn.Module):
    def __init__(
        self,
        width=64,
        depth=3,
        act=nn.Tanh(),
        quadrature_points=32,
        eps=1e-6,
    ):
        super().__init__()
        layers = []
        dims = [1] + [width] * depth + [1]
        for a, b in zip(dims[:-1], dims[1:]):
            layers.append(nn.Linear(a, b))
            if b != 1:
                layers.append(act)
        self.net = nn.Sequential(*layers)
        self.eps = float(eps)

        nodes, weights = np.polynomial.legendre.leggauss(quadrature_points)
        nodes = 0.5 * (nodes + 1.0)
        weights = 0.5 * weights
        self.register_buffer(
            "q_nodes",
            torch.tensor(nodes, dtype=torch.float32).view(-1, 1),
        )
        self.register_buffer(
            "q_weights",
            torch.tensor(weights, dtype=torch.float32).view(-1, 1),
        )

    def forward_density(self, x: torch.Tensor) -> torch.Tensor:
        return F.softplus(self.net(x)) + self.eps

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.dim() == 1:
            x = x.unsqueeze(1)
        elif x.dim() != 2 or x.shape[1] != 1:
            raise ValueError("Monotonic1D expects input with shape (N, 1) or (N,).")

        q_nodes = self.q_nodes.to(device=x.device, dtype=x.dtype)
        q_weights = self.q_weights.to(device=x.device, dtype=x.dtype)

        density_quad = self.forward_density(q_nodes)
        normalization = torch.sum(density_quad * q_weights).clamp_min(self.eps)

        mapped = x @ q_nodes.t()
        density_mapped = self.forward_density(mapped.reshape(-1, 1))
        density_mapped = density_mapped.view(x.shape[0], -1)
        integral = x * torch.sum(density_mapped * q_weights.t(), dim=1, keepdim=True)

        return integral / normalization


class DecoupledMonotonicModel(nn.Module):
    def __init__(
        self,
        width=64,
        depth=3,
        act=nn.Tanh(),
        quadrature_points=64,
        eps=1e-6,
    ):
        super().__init__()
        self.net_xi = Monotonic1D(
            width=width,
            depth=depth,
            act=act,
            quadrature_points=quadrature_points,
            eps=eps,
        )
        self.net_eta = Monotonic1D(
            width=width,
            depth=depth,
            act=act,
            quadrature_points=quadrature_points,
            eps=eps,
        )

    def forward(self, xy: torch.Tensor) -> torch.Tensor:
        x = xy[:, 0:1]
        y = xy[:, 1:2]
        xi = self.net_xi(x)
        eta = self.net_eta(y)
        return torch.cat((xi, eta), dim=1)


def _init_mlp_weights(mlp: nn.Module):
    for m in mlp.modules():
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
    last = [m for m in mlp.modules() if isinstance(m, nn.Linear)][-1]
    nn.init.zeros_(last.weight)
    nn.init.zeros_(last.bias)


def init_decoupled_model_weights(net: DecoupledMonotonicModel):
    _init_mlp_weights(net.net_xi.net)
    _init_mlp_weights(net.net_eta.net)
