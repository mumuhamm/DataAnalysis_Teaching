import matplotlib.pyplot as plt
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

def test( data):
   signal, background = data
   signal = data - background
   

s = np.arange(1,51,1)
b = 10
z1 = np.sqrt(2*((s+b)*np.log(1+s/b)-s))
z2 = s/np.sqrt(b)
z = z1/z2
sovb= s/b
#print("sqrt signal significance : ", z1, z2, z, sovb)
plt.plot(sovb, z1, color='r', label=r'$z_1 =\sqrt{2[(S+B)\ln (1+\frac{S}{B})-S]}$')
plt.plot(sovb, z2, color='g', label=r'$z+2 =\frac{S}{\sqrt{B}}$')
plt.plot(sovb, z, color='b', label=r'$z=\frac{z_1}{z_2}$')
plt.xlabel(r'$\frac{S}{B}$')
plt.ylabel("Signal significance")
plt.title("Signal signifance as a function of signal over background")
plt.legend()
  
# To load the display window
plt.show()


"""

def drumhead_height(n, k, distance, angle, t):
   kth_zero = special.jn_zeros(n, k)[-1]
   return np.cos(t) * np.cos(n*angle) * special.jn(n, distance*kth_zero)
    
    
    
theta = np.r_[0:2*np.pi:50j]
radius = np.r_[0:1:50j]
x = np.array([r * np.cos(theta) for r in radius])
y = np.array([r * np.sin(theta) for r in radius])
z = np.array([drumhead_height(1, 1, r, theta, 0.5) for r in radius])
fig = plt.figure()
ax = fig.add_axes(rect=(0, 0.05, 0.95, 0.95), projection='3d')
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_xticks(np.arange(-1, 1.1, 0.5))
ax.set_yticks(np.arange(-1, 1.1, 0.5))
ax.set_zlabel('Z')
plt.show()


# We first arrange data into a histogram
hist,bins,patches=plt.hist(obs,bins=10)

#We need to know which elements went into each bin
#np.digitize returns an array of indices corresponding to which bin from
# bins that each value in obs belongs to.
# E.g. the first entry is 5; meaning that first entry in obs belongs to
# bin number 5
bin_indices=np.digitize(obs,bins)

#under the null hypothesis, the expected values for x observations
Ex=[]
for i in range(len(bins)):
    j=np.where(bin_indices==i)[0]
    #returns the indices of obs.data which goes in bin number i
    data_in_bin=coinc[j]
    #returns the coincidence data in that bin, e.g. obs[0] is in bin 5.
    if data_in_bin.size >0:
        Ex.append(poisson.pmf(data_in_bin,MLE).sum() * sum(hist))

chisq = sum(((hist-Ex)**2)/Ex)    #this tallies up to 128.60

#   categories    params   constant
df=(len(bins)-1)   - 1     - 1

#Chi-squared statistic
p=chi2.pdf(chisq,df)     #this is 1.21 e-55
print(p)
#"automatic version"
chisquare(hist,Ex,df)   #returns p-value 2.61 e-62


def qMu(data, mu_H0, mu_H1):
  
  signal, background = data
  signalhat = data - background
  H0_data = stats.poisson.rvs(mu_H0, 1, size=data)

  lH0
  lH1=
  q0 = -2*np.log(lH0-lH1)
  retrun q0 # test statistics for discovery of a positive signal
  
  

###########################
###########################
def getPValue(data, mu_H0, mu_H1):
  apex = np.array(data)
  size = int(len(apex))
  llH0 = np.sum(mu_H0*np.log(data)+(1-mu_H0)*np.log(np.abs(1-data)))
  llH1 =  np.sum(mu_H1*np.log(data)+(1-mu_H1)*np.log(np.abs(1-data)))
  print("maximum and minimum log likelihood : ", llH0, llH1)
  qmu = 2.0*(llH1-llH0)
  print("qmu : ", qmu)
  return chi2.sf(qmu,1)
    n = 10000 # 10k simulation of the pencil length
count = 0
fraction =0

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h
sim_mean =[]
lower = []
upper = []
mean, sigma = 15, np.sqrt(12)
print(" the mean and the sigma given:  ", mean, sigma)
ci=0.607
data = []
for i in range(1,n):
  simulation = np.random.normal(mean,sigma, size=100)
  data.append(np.average(simulation))
  resampling_mean_resolution = abs(mean - np.mean(simulation))
  resampling_varinace_resolution = abs(sigma - np.std(simulation, ddof=1))
  if(resampling_mean_resolution - np.std(simulation) < 0) and (resampling_mean_resolution + np.std(simulation, ddof=1) > 0):
    fraction +=1
  xMin = simulation.min()
  xMax = simulation.max()
  afinity_outcome = np.logical_and(mean>xMin, mean<xMax)
  #print("Resampled mean and variance resolution : ", resampling_mean_resolution, resampling_varinace_resolution)
  #print(" The Maximum and the minimum value of the from the resampling", xMax, xMin)
  #print("Afinity in the outcome : ", afinity_outcome)
  m, ml, mu = mean_confidence_interval(simulation, ci)
  sim_mean.append(m)
  lower.append(ml)
  upper.append(mu)
  #print("the mean simulation, lower and upper confidence interval: ", sim_mean, lower, upper)
  resampling_standard_error = stats.sem(simulation)
  if(resampling_mean_resolution - resampling_standard_error < 0) and (resampling_mean_resolution + resampling_standard_error > 0):
    count +=1

count /=float(n)
fraction /=float(n)
print(colored("3 sigma : Fraction of regions covering the true value is?:","red"),fraction )
print(colored("1 sigma : Fraction of regions covering the true value is?:","blue"),count )
num_bins = 100
n, bins, patches = plt.hist(data[0:], num_bins, density = 1,  color ='green', alpha = 0.7)
plt.xlabel('Pencil measurement')
plt.ylabel('Frequency')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.special import gammaln # x! = Gamma(x+1)

meanlife = 550e-6
decay_lifetimes = 1/np.random.poisson((1/meanlife), size=100000)

def transformation_and_jacobian(x):
    return 1./x, 1./x**2.

def tfm_normal_pdf(x, lam):
    y, J = transformation_and_jacobian(x)
    return norm.pdf(y, lam, np.sqrt(lam)) * J

def tfm_poisson_pdf(x, mu):
    y, J = transformation_and_jacobian(x)
    # For numerical stability, compute exp(log(f(x)))
    return np.exp(y * np.log(mu) - mu - gammaln(y + 1.)) * J

hist, bins = np.histogram(decay_lifetimes, bins=50, density=True)
width = 0.8*(bins[1]-bins[0])
center = (bins[:-1]+bins[1:])/2
plt.bar(center, hist, align='center', width=width, label = 'Normalised data')

# Important: Choose a reasonable starting point
p0 = 1 / np.mean(decay_lifetimes)

norm_opt, _ = curve_fit(tfm_normal_pdf, center, hist, p0=p0)
pois_opt, _ = curve_fit(tfm_poisson_pdf, center, hist, p0=p0)

plt.plot(center, tfm_normal_pdf(center, *norm_opt), 'g--', label='(Transformed) Normal fit')
plt.plot(center, tfm_poisson_pdf(center, *pois_opt), 'r--', label='(Transformed) Poisson fit')
plt.legend(loc = 'best')
plt.tight_layout()

"""
