import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


sigma = np.sqrt(12)
print("this is the sigma vaue:  ", sigma)

def f(x1, A, B): # this is your 'straight line' y=f(x)
    return A*x1 + B
plt.figure()
m = 1.5
c = 0.5
x = np.linspace(2, 42, 20 )
y = m*x +c
print(x,y)
plt.plot(x,y,color='black', linestyle='--', linewidth=2, marker='o')
popt, pcov = curve_fit(f, x, y) # your data x, f to fit
#plt.plot(x, f(x, *popt), 'r-', label='fit: A=%5.3f, B=%5.3f' % tuple(popt))
plt.xlabel('Value of x ')
plt.ylabel('value of y = f(x) =mx+c')
plt.title('Plot of the staright line')
plt.grid(True)
plt.savefig("sl.pdf")
plt.legend()
plt.show()
  
'''m, s =  np.mean(data), np.std(data)
print("the means and sigma from the resampling: ", m, s)
s = np.random.normal(m, s, 10000)
count, bins, ignored = plt.hist(s, 20)
plt.xlabel('Pencil length')
plt.ylabel('Normalized sample frequency')
plt.title('Resampled pencil length')
plt.plot(bins, 1/(s * np.sqrt(2 * np.pi)) * np.exp( - (bins - m)**2 / (2 * s**2) ), linewidth=2, color='r')
x = np.linspace(bins.min(), bins.max(), num=100)
max_density = stats.norm.pdf(x, m, s)
scale = count.max() / max_density
plt.plot(x, scale * stats.norm.pdf(x, m, s))
'''
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

mean, sigma = 15, np.sqrt(12)
print(" the mean and the sigma given:  ", mean, sigma)
ci=0.607
for i in range(1,n):
  simulation = np.random.normal(mean,sigma, size=100)
  data.append(np.average(simulation))
  resampling_mean_resolution = abs(mean - np.mean(simulation))
  resampling_varinace_resolution = abs(sigma - np.std(simulation, ddof=1))
  if(resampling_mean_resolution - np.std(simulation, ddof=1) < 0) and (resampling_mean_resolution + np.std(simulation, ddof=1) > 0):
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
  if(resampling_mean_resolution - 2*resampling_standard_error < 0) and (resampling_mean_resolution + 2*resampling_standard_error > 0):
    count +=1
  
from random import choice
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

mu, sigma = 40, 1
population = np.random.normal(mu, sigma, 100000)

def get_confidence_interval(variance, sample_mean, sample_size, significance_level):
    margin_of_error = norm.ppf(1 - significance_level/2)*variance/np.sqrt(sample_size)
    return sample_mean - margin_of_error, sample_mean + margin_of_error

sample_size = 100
sample = np.random.choice(population, sample_size)
sample_mean = np.average(sample)
confidence_interval = get_confidence_interval(sigma, sample_mean, sample_size, 0.05)

xs = np.arange(36, 44, 0.01)
ys = norm.pdf(xs, mu, sigma)
plt.plot(xs, ys, label='population distribution')
plt.plot(confidence_interval, [0.35, 0.35], label='confidence interval for mean')
plt.legend()
plt.show()
def CI_interval(variance, sample_mean, sample_size, significance_level):
    margin_of_error = norm.ppf(1 - significance_level/2)*variance/np.sqrt(sample_size)
    return sample_mean - margin_of_error, sample_mean + margin_of_error

mean, sigma = 15, np.sqrt(12)
data_population = np.random.normal(mu, sigma, n)
sample_size = 50
sample = np.random.choice(data_population, sample_size)
sample_mean = np.mean(data_population)
sampling_mean_resolution = abs(mean - np.mean(simulation))

ci_interval = CI_interval(sigma, sample_mean, sample_size, 0.607)
xMin = data_population.min()
xMax = data_population.max()
afinity_outcome = np.logical_and(mean>xMin, mean<xMax)
print("Afinity in the outcome : ", afinity_outcome)
xs = np.arange(xMin, xMax, 0.5)
ys = norm.pdf(xs, mu, sigma)
plt.plot(xs, ys, 'o', label='population data distribution')
plt.plot(ci_interval,[0.09,0.09], label='confidence interval for mean')
plt.legend()
plt.show()

count /=float(n)
fraction /=float(n)
print(colored("3sigma : Fraction of regions covering the true value is?:","red"), count)
print(colored("1sigma : Fraction of regions covering the true value is?:","blue"), fraction)
plt.hist(simulation)

def qMu_beta(data, mu_H0, mu_H1):
  beta_base = gamma(mu_H0+mu_H1)/(gamma(mu_H0)+gamma(mu_H1))
  return beta_base*data**(mu_H0-1)*(1-data)**(mu_H1-1)
  
  N = 10000
a = np.random.normal(0, 1, N)
mean, sigma = a.mean(), a.std(ddof=1)
conf_int_a = stats.norm.interval(0.68, loc=mean, scale=sigma)

print('{:0.2%} of the single draws are in conf_int_a'
      .format(((a >= conf_int_a[0]) & (a < conf_int_a[1])).sum() / float(N)))

M = 1000
b = np.random.normal(0, 1, (N, M)).mean(axis=1)
conf_int_b = stats.norm.interval(0.68, loc=0, scale=1 / np.sqrt(M))
print('{:0.2%} of the means are in conf_int_b'
      .format(((b >= conf_int_b[0]) & (b < conf_int_b[1])).sum() / float(N)))


higgs_mass =125.35
unknown_mass = 137.25
#sample_Higgs = np.random.normal(higgs_,sigma_H0, 200)
sample_Unkwon = stats.norm.rvs(unknown_mass,sigma_H1, 100)
#data_cobined = np.concatenate(sample_Higgs, sample_Unkwon)
print("the qMu_beta case : ", qMu(sample_Unkwon, higgs_mass, unknown_mass) )


from scipy.stats import poisson
from scipy.stats import chisquare
from scipy.stats import chi2

MLE = np.mean(obs)

#H0: The data is Poisson distributed with rate lambda=MLE
#H1: The data is not Poisson distribtued

#under the null hypothesis, the expected values for x observations
# with x = 981, 982,..., 1153           (min and max of observations)

Probs = poisson.pmf(obs,MLE)
Ex = Probs * sum(obs)

chisq = sum(((obs-Ex)**2)/Ex)    #this tallies up to 253978.198

#         categories         params   constant
df=(len(np.unique(coinc)))   - 1     - 1

#Chi-squared statistic
chi2.pdf(chisq,df)     #this is 0.0

#"automatic version"
chisquare(obs,Ex,df)   #returns p-value 0.0
===================================
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import poisson

obs = np.array([1125, 1117, 1056, 1069, 1060, 1009, 1065, 1031, 1082, 1034,  985,
       1022, 1020, 1108, 1084, 1049, 1032, 1064, 1036, 1034, 1046, 1086,
       1098, 1054, 1032, 1101, 1044, 1035, 1018, 1107, 1039, 1038, 1045,
       1063,  989, 1038, 1038, 1048, 1040, 1050, 1046, 1073, 1025, 1094,
       1007, 1090, 1100, 1051, 1086, 1051, 1106, 1069, 1044, 1003, 1075,
       1061, 1094, 1052,  981, 1022, 1042, 1057, 1028, 1023, 1046, 1009,
       1097, 1081, 1147, 1045, 1043, 1052, 1065, 1068, 1153, 1056, 1145,
       1073, 1042, 1081, 1046, 1042, 1048, 1114, 1102, 1092, 1006, 1056,
       1039, 1036, 1039, 1041, 1027, 1042, 1057, 1052, 1058, 1071, 1029,
        994, 1025, 1051, 1095, 1072, 1054, 1054, 1029, 1026, 1061, 1153,
       1046, 1076])

# Simple ecdf function
def ecdf(x):
    sx = np.sort(x)
    n = sx.size
    sy = np.arange(1,n+1)/n
    return sx,sy

fig, ax = plt.subplots(figsize=(6,4))
# CDF for observed data
ecdf_x,ecdf_y = ecdf(obs)
ax.step(ecdf_x,ecdf_y,color='red',label='ECDF',
        linewidth=3,zorder=3)

# CDF for hypothetical poisson
pcdf_x = np.arange(obs.min(),obs.max()+1)
pcdf_y = 1 - poisson.cdf(obs.mean(),pcdf_x)
ax.step(pcdf_x,pcdf_y, color='k',linewidth=3,
        label=f'Poisson {obs.mean():.1f}',zorder=2)

# Random variates of same size as obs
for i in range(10):
    randp = poisson.rvs(obs.mean(),size=len(obs))
    rcdf_x,rcdf_y = ecdf(randp)
    if i == 0:
        ax.step(rcdf_x,rcdf_y, color='grey',
                label=f'Simulated Poisson',zorder=1)
    else:
        ax.step(rcdf_x,rcdf_y, color='grey',alpha=0.35,zorder=3)

ax.legend(loc='upper left')
plt.show()

================
def lill_poisson(x,sim=10000,seed=10):
    n = len(x)
    nu = np.arange(1.0,n+1)/n
    nm = np.arange(0.0,n)/n
    # Fit parameters
    m = x.mean()
    fp = poisson(m) # frozen Poisson
    # in function for KS stat
    def ks(obs):
        x = np.sort(obs)
        cv = fp.cdf(x)
        Dp = (nu - cv).max()
        Dm = (cv - nm).max()
        return np.max([Dp,Dm])
    # KS stat observed
    ks_obs = ks(x)
    # Generate simulation
    np.random.seed(seed)
    sa = np.zeros(sim)
    for i in range(sim):
        s = fp.rvs(n)
        sa[i] = ks(s)
    # calculate p-value
    p_val = np.append(sa,ks_obs).argsort()[-1]/sim
    return ks_obs, p_val, sa

kstat, pval, svals = lill_poisson(obs)
print(f'KS Stat: {kstat:0.2f}, p-value {pval:0.2f}')
==================
# Check out p-values for actual distributions
# should be approximately uniform
# On my machine, takes about 3 minutes to
# crunch through 100 simulations
# So this took me ~5 hours

res_p = []
for i in range(10000):
    if (i % 100) == 0:
        print(f'Sim {i} @ {datetime.now()}')
    mu = np.random.uniform(low=5,high=1000,size=1)[0]
    si = np.random.randint(low=10,high=100,size=1)[0]
    simloc = poisson.rvs(mu,size=si)
    ks,pv,sv = lill_poisson(simloc)
    res_p.append(pv)
# Chi-square approach bins
# Use quintiles of the hypothetical Poisson
df = pd.DataFrame(obs,columns=['counts'])
df['quantile'] = poisson.cdf(obs.mean(),obs)
df['quin'] = np.floor(df['quantile']/0.2)

obs_counts = df['quin'].value_counts()
exp_counts = len(obs)/5
chi_stat = ((obs_counts - exp_counts)**2/exp_counts).sum()
# Chi-square value of 6.66

=================




n_samples, n_features = X.shape
weights = np.zeros(n_features)

def forward(X, weights):
    return np.exp(np.dot(X, weights))

def gradient(y, y_pred, X):
    gradient = np.dot(X.T,(y - y_pred))
    return gradient


lr = 0.01
n_iter = 3000
for i in range(n_iter):
    # predicting
    y_pred = forward(X, weights)
    # computing loss
    loss = np.sqrt(np.mean((y-y_pred)**2))
    if i % 250 == 0:
        print(loss)
    
    # calculating gradient and updating weights
    calc_gradient = gradient(y, y_pred, X)
    weights -= lr * calc_gradient


import numpy as np
from scipy.stats import chisqprob

L1 = 467400. # log(likelihood) of my 1st fit
L2 = 467414. # log(likelihood) of my 2nd fit

LR = -2. * np.log(L2 / L1) # LR = -5.9905e-05

p = chisqprob(LR, 1) # L2 has 1 DoF more than L1

print 'p: %.30f' % p # p = 1.000000000000000000000000000000

five_sigma = 1 - scipy.special.erf(5 / np.sqrt(2.))                  :-)
print '5 sigma: %.30f' % five_sigma

five_sigma_check = 1 - 0.999999426696856                             :-(
print 'Check  : %.30f' % five_sigma_check


from scipy.stats.distributions import chi2
def likelihood_ratio(llmin, llmax):
    return(2*(llmax-llmin))


LR = likelihood_ratio(L1,L2)


p = chi2.sf(LR, 1) # L2 has 1 DoF more than L1

print 'p: %.30f' % p

# p: 0.000000121315450836607258011741

from scipy.stats import chi2
ll_0, ll_1 = 467400, 467414 # given, the log-likelihoods of the nested models m_0, m_1
# log likelihood for m_0 (H_0) must be <= log likelihood of m_1 (H_1)
Λ = -2 * (ll_0 - ll_1)
print(Λ)
# 28.0
df = 1 # given the difference in dof
# compute the p-value
pvalue = 1 - chi2(df).cdf(Λ) # since Λ follows χ2
print(pvalue)
# 1.2131545083660726e-07

α, df = 0.05, 1
x = np.linspace(0, 30, 1000)
plt.plot(x, chi2(df).pdf(x), label='χ2')
plt.axvline(chi2(df).ppf(1-α), color='red', label='α=0.05')
plt.scatter(Λ, 0, color='green', s=50, label='Λ')
plt.legend()
plt.title('χ2({}) distribution'.format(df), size=20)


