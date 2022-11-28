import numpy as np
from scipy.optimize import curve_fit
import math
import numpy as np

#HEP specific tools
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
from scipy import special



#Functions manipulation
from functools import partial

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
plt.rcParams.update(params)

def the_dist( mean,sigma, x):
    gaussian_dist = (1/np.sqrt(2*np.pi*(sigma**2)))*np.exp(-(mean-x)**2/2 * sigma**2)
    return gaussian_dist 

data_x = stats.poisson.rvs(20, 6, 10000)
mean = 5.6
sigma = 0.15
print("Gaussian points : ", the_dist(mean, sigma, data_x)) 
