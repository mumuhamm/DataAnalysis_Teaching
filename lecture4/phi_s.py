#!/usr/bin/env python3
"""
Compare confidence intervals with cuts
======================================

Using the b-physics phi_s variable, show difference in coverage between
intersection of 5 confidence intervals and a combined confidence interval.
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from phi_s_zoom import zoom_effect
from style import make_style


# True value of \phi_s - use SM prediction
# https://arxiv.org/pdf/1909.12524.pdf
# eq. 91
truth = -0.0369

# \phi_s measurement standard deviations
# https://hflav-eos.web.cern.ch/hflav-eos/osc/summer_2017/HFLAV_phis_inputs.pdf
sigma_data = {"atlas_combined": 0.078,
              "cms": 0.097,
              "lhcb_combined": 0.037,
              "lhcb_psi_2S_phi": 0.285,  # averaged asymmetric +0.29/-0.28
              "lhcb_DD": 0.17}
sigma = np.array(list(sigma_data.values()))

# Weighted least squares parameters
weights = 1. / sigma**2
combined_sigma = weights.sum()**-0.5
normalized_weights = weights / weights.sum()

# Confidence level etc
alpha = 0.05
z = norm.isf(alpha * 0.5)
beta = 1. - alpha
beta_percent = 100. * beta

# Number of repeats etc
n_measurements = len(sigma)
n_trials = 10000
n_show = 100
seed = 125


def confidence_interval(data):
    """
    Unknown mean, known variance. Weighted least squares.

    @returns Confidence interval
    """
    sample_mean = np.dot(normalized_weights, data)
    return (sample_mean - z * combined_sigma,
            sample_mean + z * combined_sigma)

def individual_confidence_intervals(data):
    """
    @returns Collection of individual intervals from each group of draws from Gaussian
    """
    u = data + z * sigma
    l = data - z * sigma
    return np.column_stack((l, u))

def cuts_interval(data):
    """
    @returns Interval from cut on each draw from Gaussian.
    """
    max_ = min(data + z * sigma)
    min_ = max(data - z * sigma)
    if max_ > min_:
        return (min_, max_)
    return (np.nan, np.nan)

def pseudo_data():
    """
    @returns Set of independent draws from Gaussian distribution
    """
    return norm.rvs(truth, sigma)

def contains_fraction(intervals):
    """
    @returns Fraction of times intervals contain true value
    """
    contains = sum(min(i) <= truth <= max(i) for i in intervals)
    return float(contains) / len(intervals)


if __name__ == "__main__":

    # Set random seed for reproducibility
    np.random.seed(seed)

    # Set plot styling
    make_style()

    # Pseudo-data for each panel
    data = [pseudo_data() for _ in range(n_trials)]

    conf = [confidence_interval(d) for d in data]
    ind_conf = [individual_confidence_intervals(d) for d in data]

    conf_beta = contains_fraction(conf)
    print("Confidence intervals contained truth in {} pseudo-experiments. Expected {}.".format(conf_beta, beta))

    cuts = [cuts_interval(d) for d in data]
    cuts_beta = contains_fraction(cuts)
    cuts_beta_theory = beta**n_measurements
    print("Cuts contained truth in {} pseudo-experiments. Expected {}.".format(cuts_beta, cuts_beta_theory))

    # Color scheme

    c_inside = '#4da6ff'
    c_individual = 'k'
    c_outside = '#ff531a'
    c_empty = '0.8'

    # Arrange axes

    fig = plt.figure()
    grid = gridspec.GridSpec(3, 4, width_ratios=[0.7, 1, 1, 0.7],
                             height_ratios=[0.1, 1, 0.1], wspace=0.3)

    ax1 = plt.subplot(grid[1, 0])
    ax2 = plt.subplot(grid[:, 1])
    ax3 = plt.subplot(grid[:, 2])
    ax4 = plt.subplot(grid[1, 3])
    ax = [ax1, ax2, ax3, ax4]

    # Show repeats in central panels

    for i, interval in enumerate(conf[:n_show]):
        c = c_inside if min(interval) <= truth <= max(interval) else c_outside
        ax[1].plot(interval, [i] * 2, color=c)

    for i, interval in enumerate(cuts[:n_show]):
        c = c_inside if min(interval) <= truth <= max(interval) else c_outside
        if np.any(np.isnan(interval)):
            # empty interval
            ax[2].axhline(i, color=c_empty)
        else:
            ax[2].plot(interval, [i] * 2, color=c)

    # Make zoom data

    zoom = (46, 48)
    indexes = np.arange(zoom[0], zoom[1] + 1, 1)
    pad_between_data = 5
    y = 0

    for index in indexes:
        for j, interval in enumerate(ind_conf[index]):
            c = c_inside if min(interval) <= truth <= max(interval) else c_outside
            center = np.mean(interval)
            error = interval[1] - center
            ax[0].errorbar(center, y + j, xerr=error, color=c_individual, capsize=1.5, lw=1.)
            ax[3].errorbar(center, y + j, xerr=error, color=c_individual, capsize=1.5, lw=1.)

        for ax_index in [0, 3]:
            interval = conf[index] if ax_index == 0 else cuts[index]
            empty = np.any(np.isnan(interval))
            c = c_empty if empty else c_inside if min(interval) <= truth <= max(interval) else c_outside
            center = truth if empty else np.mean(interval)
            error = 2. if empty else interval[1] - center
            ax[ax_index].errorbar(center, y + 0.5 * len(indexes) + 0.5,
                                  xerr=error, color=c, capsize=24, capthick=1.5, lw=40., zorder=-1)

        y += len(indexes) + pad_between_data

    # Set identical limits on central plots

    l = min(ax[1].get_xlim()[0], ax[2].get_xlim()[0])
    u = max(ax[1].get_xlim()[1], ax[2].get_xlim()[1])

    for a in ax[1:3]:
        a.set_ylim(-1., n_show)

    for a in ax:
        a.axvline(truth, color='k', lw=1.5, label="Truth (SM prediction)")
        a.set_xlim(l, u)
        a.set_yticks([])

    for a in (ax[0], ax[3]):
        a.set_xticks([])
        a.set_xlim(3 * l, 3 * u)
        a.set_ylim(-0.5 * pad_between_data, y - 0.5 * pad_between_data)

    # Instead of set_title to center it on two axes
    ax[1].text(ax[1].get_xlim()[1], ax[1].get_ylim()[1] + 2.,
               "Combined ${0:.0f}\%$ confidence intervals \n \\textit{{Contains true value: ~~~~~~~~~~~~~~~~~{1:.0f}\%}}"
               .format(100. * beta, 100. * conf_beta),
               horizontalalignment='right', fontsize=12)
    ax[2].text(ax[2].get_xlim()[0], ax[2].get_ylim()[1] + 2.,
               "Intersection of ${0:.0f}\%$ confidence intervals \n \\textit{{{1:.0f}\%}}"
               .format(100. * beta, 100. * cuts_beta),
               horizontalalignment='left', fontsize=12)

    # Instead of set_xlabel to center it on two axes
    ax[1].text(0., -15., r"Phase $\phi_s$ (radians)")

    # Annotations etc
    ax[0].text(-1., 5., r"Five measurements", rotation=90, fontsize=11)
    ax[0].annotate(r"$\{$", fontsize=40,
                   xy=(0.05, 0.505), xycoords='figure fraction')
    ax[1].text(0.155, 37, r"100 repetitions", rotation=90, fontsize=11)

    zoom_effect(fig, ax, zoom[0] - 0.5, zoom[1] + 0.5)

    plt.savefig("phi_s.pdf", bbox_inches='tight')
