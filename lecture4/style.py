#!/usr/bin/env python3
"""
Common styling for all figures
==============================
"""

from matplotlib import rc


def make_style():
    """
    Applies style commands to global state
    """
    rc('text', usetex=True)
    rc('text.latex', preamble=(r'\setlength\parindent{0pt}'
                               r'\usepackage[scaled=0.92]{helvet}'
                               r'\usepackage[helvet]{sfmath}'
                               r'\usepackage{bm}'
                               r'\usepackage{siunitx}'
                               r'\sisetup{group-separator={,},group-minimum-digits={3}}'))
    rc('font', **{'family': 'sans-serif', 'size': 14})
    rc('axes', **{'grid': False, 'titlesize': 14, 'labelsize': 14})
    rc('figure', **{'titlesize': 14, 'dpi': 600})
    rc('legend', **{'fontsize': 12, 'title_fontsize': 14, 'handlelength': 1., 'frameon': False})
