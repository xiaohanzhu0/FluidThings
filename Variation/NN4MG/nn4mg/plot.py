import matplotlib.pyplot as plt


def plot_grid(X, Y, title=None):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.plot(X, Y, "k", linewidth=0.1)
    ax.plot(X.T, Y.T, "k", linewidth=0.1)
    ax.set_aspect("equal", adjustable="box")
    if title:
        ax.set_title(title)
    return fig, ax
