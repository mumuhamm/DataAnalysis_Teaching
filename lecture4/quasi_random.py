#!/usr/bin/env python3
"""
Benefit of quasi-random over random and grid search
===================================================
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc

from qmcpy import Sobol
from style import make_style


ln2_npoints = 8
npoints = 2**ln2_npoints
np.random.seed(125)


def peak(x, mu, sigma):
    """
    @returns Gaussian peak
    """
    return np.exp(-(x - mu)**2 / (2. * sigma**2))

def important(x):
    """
    @brief Function of important variable; two peaks
    """
    return peak(x, 0.7, 0.3) + 1.4 * peak(x, 0.3, 0.1)

if __name__ == "__main__":

    # Set plot styling
    make_style()
    rc("lines", markersize=5)
    rc('font', **{'size': 22})
    rc('axes', **{'titlesize': 22, 'labelsize': 22})
    rc('figure', **{'titlesize': 22})

    # Create figure
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(14, 7), gridspec_kw={'height_ratios': [1, 3]})

    x = np.linspace(0, 1., 100)
    y = important(x)
    for a in ax[0]:
        a.plot(x, y)

    # grid
    points_per_dim = int(2**(ln2_npoints / 2))
    points = [points_per_dim] * 2
    grid = np.array(list(np.ndindex(*points)))
    grid = grid / grid.max()
    ax[1, 0].scatter(grid[:, 0], grid[:, 1])
    ax[0, 0].set_title("Grid")
    ax[0, 0].scatter(grid[:, 0], important(grid[:, 0]))

    # random
    sample = np.random.rand(npoints, 2)
    ax[1, 1].scatter(sample[:, 0], sample[:, 1])
    ax[0, 1].set_title("Random")
    ax[0, 1].scatter(sample[:, 0], important(sample[:, 0]))

    # quasi-random - sobol
    sampler = Sobol(dimension=2, randomize=True, seed=7)
    qsample = sampler.gen_samples(npoints)
    ax[1, 2].scatter(qsample[:, 0], qsample[:, 1])
    ax[0, 2].set_title("Quasi-random")
    ax[0, 2].scatter(qsample[:, 0], important(qsample[:, 0]))

    for a in ax[1]:
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        a.get_xaxis().set_ticks([])
        a.get_yaxis().set_ticks([])

    ax[1, 1].set_xlabel('Important parameter', labelpad=30)
    ax[1, 0].set_ylabel('Un-important parameter', labelpad=30)
    ax[0, 0].set_ylabel('Function', labelpad=30)

    for a in ax[0]:
        for side in ['top', 'right', 'bottom', 'left']:
            a.spines[side].set_visible(False)
            a.tick_params(axis='both', which='both', labelbottom=False, labelleft=False, bottom=False, left=False)

    plt.tight_layout()
    plt.savefig("quasi_random.pdf")
