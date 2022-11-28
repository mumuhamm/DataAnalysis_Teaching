#!/usr/bin/env python3
"""
Zoom effect
===========

Functions to create a zoom effect for the phi_s figure
"""

from matplotlib.patches import Rectangle, ConnectionPatch


def zoom_effect(fig, ax, ymin, ymax, color="teal", alpha=0.1):
    """
    Adds shading and lines to indicate zoomed regions
    """
    # Shade entire side panels

    for a in (ax[0], ax[3]):
        a.patch.set_facecolor(color)
        a.spines['bottom'].set_color(color)
        a.spines['top'].set_color(color)

        a.patch.set_alpha(alpha)
        a.spines['bottom'].set_alpha(0.)
        a.spines['top'].set_alpha(0.)

    # Shade patch of central panels

    for a in (ax[1], ax[2]):
        r = Rectangle((a.viewLim.x0, ymin), a.viewLim.x1 - a.viewLim.x0, ymax - ymin,
                      color=color, lw=0, alpha=alpha)
        a.add_patch(r)

    # Draw lines from central panels to zoomed regions

    add_line(fig, ax[1], ax[0], ax[1].viewLim.x0, ymin, ax[0].viewLim.x1, ax[0].viewLim.y0)
    add_line(fig, ax[1], ax[0], ax[1].viewLim.x0, ymax, ax[0].viewLim.x1, ax[0].viewLim.y1)
    add_line(fig, ax[2], ax[3], ax[2].viewLim.x1, ymin, ax[3].viewLim.x0, ax[3].viewLim.y0)
    add_line(fig, ax[2], ax[3], ax[2].viewLim.x1, ymax, ax[3].viewLim.x0, ax[3].viewLim.y1)

def add_line(fig, ax1, ax2, x1, y1, x2, y2):
    """
    Adds line across figure
    """
    trans1 = ax1.transData
    trans2 = ax2.transData
    inv_trans = fig.transFigure.inverted()

    a = inv_trans.transform(trans1.transform((x1, y1)))
    b = inv_trans.transform(trans2.transform((x2, y2)))

    patch = ConnectionPatch(xyA=a, coordsA=fig.transFigure, xyB=b,
                            coordsB=fig.transFigure, linestyle="--")
    fig.patches.append(patch)
