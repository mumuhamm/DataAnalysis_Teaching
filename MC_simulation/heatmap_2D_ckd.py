import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial import cKDTree


def data_coord2view_coord(p, resolution, pmin, pmax):
    dp = pmax - pmin
    dv = (p - pmin) / dp * resolution
    return dv


n = 1000
xs = np.random.randn(n)
ys = np.random.randn(n)

resolution = 250

extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]
xv = data_coord2view_coord(xs, resolution, extent[0], extent[1])
yv = data_coord2view_coord(ys, resolution, extent[2], extent[3])


def kNN2DDens(xv, yv, resolution, neighbours, dim=2):
    """
    """
    # Create the tree
    tree = cKDTree(np.array([xv, yv]).T)
    # Find the closest nnmax-1 neighbors (first entry is the point itself)
    grid = np.mgrid[0:resolution, 0:resolution].T.reshape(resolution**2, dim)
    dists = tree.query(grid, neighbours)
    # Inverse of the sum of distances to each grid point.
    inv_sum_dists = 1. / dists[0].sum(1)

    # Reshape
    im = inv_sum_dists.reshape(resolution, resolution)
    return im


fig, axes = plt.subplots(2, 2, figsize=(8, 8))
for ax, neighbours in zip(axes.flatten(), [0, 16, 32, 63]):

    if neighbours == 0:
        ax.plot(xs, ys, 'k.', markersize=5)
        ax.set_aspect('equal')
        ax.set_title("Scatter Plot")
    else:

        im = kNN2DDens(xv, yv, resolution, neighbours)

        ax.imshow(im, origin='lower', extent=extent, cmap=cm.Blues)
        ax.set_title("Smoothing over %d neighbours" % neighbours)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

plt.savefig('new.png', dpi=150, bbox_inches='tight')
plt.show()
