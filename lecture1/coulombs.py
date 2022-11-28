import matplotlib.pyplot as plt
import numpy as np
from math import *
from copy import *


def efield(xs, ys):
    qpos = {1: (-1,0), -1: (1,0)}
    n = len(xs)
    exs = [[0. for k in range(n)]for j in range(n)]
    eys = deepcopy(exs)
    for j, x in enumerate(xs):
        for k, y in enumerate(ys):
            for q, pos in qpos.items():
                posx , posy = pos
                r = sqrt((x-posx)**2 + (y-posy)**2)
                exs[k][j] += q*(x-posx)/r**3
                eys[k][j] += q*(y-posy)/r**3
    return exs, eys

def plotefield(boxl, n):
    xs = [-boxl+i*2*boxl/(n-1) for i in range(n)]
    ys = xs[:]
    exs, eys = efield(xs, ys)
    xs = np.array(xs)
    ys = np.array(ys)
    exs = np.array(exs)
    eys = np.array(eys)
    plt.streamplot(xs, ys, exs, eys, density=1.5, color='m')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()

plotefield(2.0, 20)
