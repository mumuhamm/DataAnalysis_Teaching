import matplotlib.pyplot as plt
from math import *
import numpy as np

def exact_f(x):
    return (exp(x)-1)/x

def trick(x):
    w = exp(x)
    if w == 0.:
        val = -1/x
    elif w == 1.0:
        val = 1.
    else:
        val = w -1 /log(w)
    return val

xs = [10**(-i) for i in ( 12, 13, 14)]
xs += [-10**(-i) for i in (11, 12, 14)]
fvals = [ exact_f(x) for x in xs]
gvals = [trick(x) for x in xs]
for x, fval , gval in zip(xs, fvals, gvals):
    print("[", x, fval, gval, "]")
