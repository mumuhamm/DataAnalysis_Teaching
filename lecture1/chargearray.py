from kahansum import kahansum as ks
from scipy import * 
from math import * 
import numpy as np 

def chargearray(nvals):
    vals = [-0.5+i/(nvals-1) for i in range(nvals)]
    qpos = {}
    for i, posx in enumerate(vals):
        for j, posy in enumerate(vals):
            count = i*nvals + j +1 
            key = 1.02*count if (i+j)%2 == 0 else -count 
            qpos[key] = posx, posy
    return qpos

def vecmag(rs):
    sq = [r**2 for r in rs]
    return np.sqrt(ks(sq))

def potential( qpos, rs):
    potvals = []
    for q, pos in qpos.items():
        diffs = [r - po for r, po in zip(rs,pos)]
        R = vecmag(diffs)
        potvals.append(q/R)
    return ks(potvals)

if __name__ == '__main__':
    for i in range(2, 40):
        qpos = chargearray(i)
        for y in 3,-3:
            rs = [0., y]
            potval = potential(qpos, rs)
            print(potval)
    
