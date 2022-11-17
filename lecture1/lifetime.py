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
plt.show()
