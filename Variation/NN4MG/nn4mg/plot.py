import matplotlib.pyplot as plt


def plot_grid(X, Y, title=None):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.plot(X, Y, "k", linewidth=0.1)
    ax.plot(X.T, Y.T, "k", linewidth=0.1)
    ax.set_aspect("equal", adjustable="box")
    if title:
        ax.set_title(title)
    return fig, ax


def plot_scalar_field(X, Y, field, title=None, cmap="viridis"):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    mesh = ax.pcolormesh(X, Y, field, shading="auto", cmap=cmap)
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(mesh, ax=ax, shrink=0.8)
    if title:
        ax.set_title(title)
    return fig, ax
