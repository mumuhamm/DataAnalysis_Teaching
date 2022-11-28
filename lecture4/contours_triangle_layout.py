#!/usr/bin/env python3
"""
Contours figure
===============

Figure showing individual 2d likelihoods, the intersection of the
individual 95% regions and a combined 95% region on a 2d plane.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import norm, chi2

from style import make_style


# Examples of 95% contours on a 2d plane

def limit_a(x):
    """
    @returns Quadratic lower limit on y
    """
    l = 2. * (0.4**2 - x**2)
    l[l < 0.] = 0.
    return l

def limit_b(y):
    """
    @returns Lower limit on x (in this case indepdent of y)
    """
    return 0.1

def limit_c(x):
    """
    @returns Cubic upper limit on y
    """
    return - ((x - 0.7) / 0.25)**3 + 0.75

def limit_d(x):
    """
    @returns Straight line lower limit on y
    """
    return x - 0.5

class Smear(object):
    """
    Functor to smear limit with a Gaussian.

    We ensure chi-squared difference corresponds to limit.
    """
    def __init__(self, limit=0., sigma=0.2, level=0.95, dof=2):
        delta_chi_squared = chi2.ppf(level, dof)
        r = np.exp(-0.5 * delta_chi_squared)
        norm_ = norm(limit, sigma)
        self.displace = norm_.ppf(r)
        self.logcdf = norm_.logcdf

    def __call__(self, x):
        """
        @returns Log-likelihood
        """
        return self.logcdf(x + self.displace)

smear = Smear()

def combined_loglike(x, y):
    """
    @returns Possible combined log-likelihood corresponding to the limits
    """
    d1 = y - limit_a(x)
    d2 = x - limit_b(y)
    d3 = limit_c(x) - y
    d4 = y - limit_d(x)
    return smear(np.array([d1, d2, d3, d4])).sum(axis=0)

def loglike_a(x, y):
    """
    @returns Possible log-likelihood corresponding to limit a
    """
    d1 = y - limit_a(x)
    return smear(np.array([d1])).sum(axis=0)

def loglike_b(x, y):
    """
    @returns Possible combined log-likelihood corresponding to limit b
    """
    d4 = y - limit_d(x)
    return smear(np.array([d4])).sum(axis=0)

def loglike_c(x, y):
    """
    @returns Possible combined log-likelihood corresponding to limit c
    """
    d2 = x - limit_b(y)
    return smear(np.array([d2])).sum(axis=0)

def loglike_d(x, y):
    """
    @returns Possible combined log-likelihood corresponding to limit d
    """
    d3 = limit_c(x) - y
    return smear(np.array([d3])).sum(axis=0)

def combined_contours(x, y):
    """
    @returns Whether allowed by crude combination of the limits
    """
    d1 = y > limit_a(x)
    d2 = x > limit_b(x)
    d3 = limit_c(x) > y
    d4 = y > limit_d(x)
    combined = np.logical_and(np.logical_and(d1, d2), np.logical_and(d3, d4))
    return combined


if __name__ == "__main__":

    # Set plot styling
    make_style()

    # Create figure
    fig = plt.figure(constrained_layout=False, figsize=(7, 7))

    # Add invisible axis instance to draw the arrows
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], frame_on=False)
    ax.set_xlim([0., 1.])
    ax.set_ylim([0., 1.])
    ax.set_yticks([])
    ax.set_xticks([])
    plt.arrow(0.2, 0.53, -0.05, -0.06, width=0.012, head_length=0.03, facecolor='0.60', edgecolor='0.05')
    plt.arrow(1-0.2, 0.53, 0.05, -0.06, width=0.012, head_length=0.03, facecolor='0.60', edgecolor='0.05')

    #
    # Top panel
    #

    panel_width = 0.30
    panel_height = panel_width
    ax = fig.add_axes([0.35, 0.58, panel_width, panel_height])

    N = 500
    x = np.linspace(0., 1., N)
    g = np.meshgrid(x, x)

    loglike = loglike_a(g[0], g[1])
    delta_chi_squared_a = -2. * (loglike - loglike.max())
    levels_a = np.arange(1.5, 12, 1.5)
    c = ax.contourf(x, x, delta_chi_squared_a, levels=levels_a, cmap=cm.Blues, extend='max', alpha=1.0)

    loglike = loglike_b(g[0], g[1])
    delta_chi_squared_b = -2. * (loglike - loglike.max())
    levels_b = np.arange(1.5, 12, 1.5)
    c = ax.contourf(x, x, delta_chi_squared_b, levels=levels_b, cmap=cm.Reds, extend='max', alpha=0.6)

    loglike = loglike_d(g[0], g[1])
    delta_chi_squared_d = -2. * (loglike - loglike.max())
    levels_d = np.arange(1.5, 12, 1.5)
    c = ax.contourf(x, x, delta_chi_squared_d, levels=levels_d, cmap=cm.Greens, extend='max', alpha=0.6)

    loglike = loglike_c(g[0], g[1])
    delta_chi_squared_c = -2. * (loglike - loglike.max())
    levels_c = np.arange(1.5, 12, 1.5)
    c = ax.contourf(x, x, delta_chi_squared_c, levels=levels_c, cmap=cm.Oranges, extend='max', alpha=0.4)

    ax.set_title(r"Individual likelihood functions", fontsize=12)
    ax.set_xlim(0., 1)
    ax.set_ylim(0., 1)
    ax.set_xlabel("Model parameter $x$")
    ax.set_ylabel("Model parameter $y$")
    ax.set_aspect(1)


    #
    # Bottom left panel
    #

    ax = fig.add_axes([0.1, 0.1, panel_width, panel_height])

    N = 1000
    x = np.linspace(0., 1., N)
    style = {'hatch': "x", 'lw': 2, 'alpha': 0.5}

    ax.fill_between(x, limit_a(x), **style)
    ax.fill_between(x, 1., where=x < limit_b(x), **style)
    ax.fill_between(x, limit_c(x), 1., **style)
    ax.fill_between(x, limit_d(x), **style)

    N = 500
    x = np.linspace(0., 1., N)
    g = np.meshgrid(x, x)
    contour = combined_contours(g[0], g[1])
    ax.contour(x, x, contour, colors="0.2", linewidths=3, linestyles="--")

    ax.set_title(r"Simplistic: overlay $95\%$ regions", fontsize=12)
    ax.set_xlim(0., 1)
    ax.set_ylim(0., 1)
    ax.set_xlabel("Model parameter $x$")
    ax.set_ylabel("Model parameter $y$")
    ax.set_aspect(1)

    #
    # Bottom right panel
    #

    ax = fig.add_axes([0.6, 0.1, panel_width, panel_height])
    cax = fig.add_axes([0.6 + 1.03 * panel_height, 0.1, 0.02, panel_height])

    N = 500
    x = np.linspace(0., 1., N)
    g = np.meshgrid(x, x)

    loglike = combined_loglike(g[0], g[1])
    delta_chi_squared = -2. * (loglike - loglike.max())
    levels = np.arange(0, 20, 2.0)
    c = ax.contourf(x, x, delta_chi_squared, levels=levels, cmap=cm.Blues, extend='max')

    cbar = fig.colorbar(c, cax=cax, ticks=np.arange(0, 20, 2.0))
    cbar.ax.set_title(r"$-2\ln\mathcal{L}$", loc='left')

    # Contour for combined log-likelihood on right panel
    ax.contour(x, x, delta_chi_squared,
               levels=[chi2.isf(0.05, 2)], colors="black", linewidths=3.5, zorder=3)
    ax.contour(x, x, delta_chi_squared,
               levels=[chi2.isf(0.05, 2)], colors="red", linewidths=3, zorder=3)

    # Contour from left panel
    contour = combined_contours(g[0], g[1])
    ax.contour(x, x, contour, colors="0.2", linewidths=3, linestyles="--")

    ax.set_title(r"Better: combine likelihoods", fontsize=12)
    ax.set_xlim(0., 1)
    ax.set_ylim(0., 1)
    ax.set_xlabel("Model parameter $x$")
    ax.set_ylabel("Model parameter $y$")
    ax.set_aspect(1)

    # Save figure
    plt.savefig("contours_triangle_layout.pdf", bbox_inches='tight')
