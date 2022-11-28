import numpy as np
from math import * 
from scipy import *

def kahansum(xs):
    s =0 ; e =0 
    for x in xs:
        temp = s
        y = x + e
        s = temp + y 
        e = (temp -s ) + y
    return s
if __name__ == '__main__':
    xs = [0.1, 0.4, 0.6]
    print(xs)
    print(np.sum(xs), kahansum(xs))
