import torch
import torch.nn as nn


class Model(nn.Module):
    """
    Residual network R_theta: (x,y) -> (Delta xi, Delta eta).
    Full map is: s(x,y) = (x,y) + g(x,y) * R_theta(x,y).
    """

    def __init__(self, width=128, depth=4, act=nn.Tanh()):
        super().__init__()
        layers = []
        in_dim, out_dim = 2, 2

        dims = [in_dim] + [width] * depth + [out_dim]
        for a, b in zip(dims[:-1], dims[1:]):
            layers.append(nn.Linear(a, b))
            if b != out_dim:
                layers.append(act)

        self.net = nn.Sequential(*layers)

    def forward(self, xy: torch.Tensor) -> torch.Tensor:
        x = xy[:, 0:1]
        y = xy[:, 1:2]
        dirichlet1 = x * (1.0 - x)
        dirichlet2 = y * (1.0 - y)
        return torch.cat(
            (
                x + dirichlet1 * self.net(xy)[:, 0:1],
                y + dirichlet2 * self.net(xy)[:, 1:2],
            ),
            dim=1,
        )


def init_model_weights(net: nn.Module):
    for m in net.modules():
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)

    last = [m for m in net.modules() if isinstance(m, nn.Linear)][-1]
    nn.init.zeros_(last.weight)
    nn.init.zeros_(last.bias)
