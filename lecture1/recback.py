from math import * 
import matplotlib.pyplot as plt
import numpy as np

def backward(nmax =31):
    oldint = 0.01 
    for n in reversed( range(20, nmax)):
        print(n, oldint)
        newint = (oldint + exp(-1))/n
        oldint = newint 
print(" n=20, answer = 0.018")
print("n, f(n)")
backward()
