#Color printing
from termcolor import colored

#General data operations library
import math
import numpy as np
import os, sys
import matplotlib as mp
import pylab as pl
#HEP specific tools
import pandas as pd
import scipy as sp
import scipy.constants as scipy_constants
from scipy.stats import poisson
from scipy.stats import chi2
from scipy.stats import norm
from scipy.stats import gamma
import scipy.stats as stats
from random import choice
from scipy.stats.distributions import chi2
#Plotting libraries
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#Functions manipulation
from functools import partial
#Increase plots font size
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
plt.rcParams.update(params)
plt.rc('legend',fontsize=10)
"""qmu = -2ln(ll(H0)/ll(H1)) test statistics musig =25, mubkg = 200,
likelihoods are the product of the poisson processes,
log likelihood is the sum of the processes
Produce test statitics using qmu, calculate p value"""

muSgn = 25
muBkg = 200
muObs = muBkg + muSgn
deltamuBkg = 10 #5% fluctuation in bkg estimation , for the sake of generality
eff_muBkg = (muBkg/deltamuBkg)**2
scale_factor = eff_muBkg/muBkg
print ("Assumed effective background count (eff_muBkg):        %10.2f" % eff_muBkg)
print ("Assumed effective background scale factor (scale_factor): %10.2f" % scale_factor)

#marginalized likelihood
def Int_likelihood(muObs, eff_muBkg, scale_factor, s):
    u = 1.0/(1 + scale_factor)
    if type(muObs) == type(0):
        r = np.arange(0, muObs+1)
        beta = stats.beta.pdf(u, r + 1, eff_muBkg)
        poisson = stats.poisson.pmf(muObs - r, s)
        prob = (1 - u) * (1 - u) * sum(beta * poisson) / eff_muBkg
        return prob
    else:
        prob = []
        for n in muObs:
            r = np.arange(0, n+1)
            beta = stats.beta.pdf(u, r + 1, eff_muBkg)
            poisson = stats.poisson.pmf(n - r, s)
            prob.append((1 - u) * (1 - u) * sum(beta * poisson) / eff_muBkg)
        return np.array(prob)
         
#profile likelihood

def Prof_likelihood(muObs, eff_muBkg, k, s):
    if type(muObs) == type(0):
        y = muObs + eff_muBkg - (1+scale_factor)*s
        b = y + np.sqrt(y*y+4*(1+scale_factor)*eff_muBkg*s)
        b /= 2*(1+scale_factor)
        prob1 = sp.stats.gamma.pdf(s + b, muObs + 1)
        prob2 = sp.stats.gamma.pdf(scale_factor*b, eff_muBkg + 1)
        return prob1 * prob2
    else:
        prob = []
        for n in muObs:
            y = n + eff_muBkg - (1+scale_factor)*s
            b = y + np.sqrt(y*y+4*(1+scale_factor)*eff_muBkg*s)
            b /= 2*(1+scale_factor)
            prob1 = sp.stats.gamma.pdf(s + b, muObs + 1)
            prob2 = sp.stats.gamma.pdf(scale_factor*b, eff_muBk + 1)
            prob.append(prob1*prob2)
        return np.array(prob)
  
#plot the likelihoods

def plotLikelihoods(Prof_likelihood, Int_likelihood, muObs, eff_muBkg, scale_factor):
    plt.figure(figsize=(6, 6))
    x  = np.arange(0, 100, 0.1) # x = [0, 0.1, ...]
    u1 = [Int_likelihood(muObs, eff_muBkg, scale_factor, s) for s in x]; u1 = u1 / max(u1)
    u2 = [Prof_likelihood(muObs, eff_muBkg, scale_factor, s)  for s in x]; u2 = u2 / max(u2)
    plt.plot(x, u1,
             color=(0,0,1),
             linewidth=1,
             label=r'$P(DATA|s)$')
    plt.fill_between(x, u1, alpha=0.05, color=(0,0,1))
    pl.legend()
    plt.plot(x, u2,
             color=(1,0,0),
             linewidth=2,
             label=r'$\mathcal{L}_p(s)=P(DATA|s, \hat{b}=g(s))$')
    pl.legend()
    plt.legend(loc='upper right', title=r'$\mathcal{L}-metric$')
    axes = plt.gca()
    umin, umax = axes.get_ylim()
    delta = 0.2
    ii = int(1.3*umax/delta)
    umax = ii * delta
    axes.set_ylim((0, umax))
    axes.set_xlabel(r'$s$', fontsize=15)
    axes.set_ylabel(r'$P(DATA|s)$', fontsize=15)
    plt.autoscale()
    plt.savefig("fig_profile_vs_marginal.png")
    plt.show()
  
plotLikelihoods(Prof_likelihood, Int_likelihood, muObs, eff_muBkg, scale_factor)

# Test statistics generation 
def q_Mu(muObs, eff_muBkg, scale_factor, muSgn, s):
    L_H0 = Prof_likelihood(muObs, eff_muBkg, scale_factor, 0)
    L_H1 = Prof_likelihood(muObs, eff_muBkg, scale_factor, muSgn)
    test_stat= -2*np.log(L_H0/L_H1)
    Z  = np.sqrt(test_stat)
    return (test_stat, Z)
print ("test_stats_obs(0), Z = %10.2f, %10.2f" % q_Mu(muObs, eff_muBkg, scale_factor, muSgn, 0))


# pvalue calculation
def pval_nouncer(muObs, muBkg):
    p_val = sp.special.gammainc(muObs, muBkg)
    Z  = sp.special.erfinv(1-2*p_val)*np.sqrt(2)
    return (p_val, Z)

def pval_uncer(muObs, eff_muBkg, scale_factor):
    n  = np.arange(0, muObs)
    p_val = 1 - sum(Int_likelihood(n, eff_muBkg, scale_factor, 0))
    Z  = sp.special.erfinv(1-2*p_val)*np.sqrt(2)
    return (p_val, Z)
    
    
    
pval_noErr = pval_nouncer(muObs, muBkg)
pval_withErr = pval_uncer(muObs, eff_muBkg, scale_factor)
print ("p-value without uncertainty ,   Z = %10.3e, %10.2f" % pval_noErr)
print ("p-value including uncertainty , Z = %10.3e, %10.2f" % pval_withErr)
