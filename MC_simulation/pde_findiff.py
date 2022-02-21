import math
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

T = 1  #Time to Expiry in Years
E = 100  #Strike
r = .05  #Risk Free Rate
SIGMA = .20  #Volatility
Type = True   #Type of Option True=Call False=Put
NAS = 40  #Number of Asset Steps - Higher is more accurate, but more time consuming

ds = 2 * E / NAS  #Asset Value Step Size
dt = (0.9/NAS/NAS/SIGMA/SIGMA)  #Time Step Size
NTS = int(T / dt) + 1  #Number of 

value_matrix = np.zeros((int(NAS+1), int(NTS)))
asset_price = np.arange(NAS*ds,-1,-ds)
value_matrix[:,-1]= np.maximum(asset_price - E,0)

for x in range(1,NTS):
    value_matrix[-1,-x-1] = value_matrix[-1,-x]* math.exp(-r*dt)
for x in range(1,int(NTS)):
    for y in range(1,int(NAS)):
        #Evaluate Option Greeks
        Delta = (value_matrix[y-1,-x] - value_matrix[y+1,-x]) / 2 / ds
        value_matrix[y+1,-x]
        Gamma = (value_matrix[y-1,-x] - (2 * value_matrix[y,-x]) + value_matrix[y+1,-x]) / ds / ds
        Theta = (-.5 * SIGMA**2 * asset_price[y]**2 * Gamma) - (r * asset_price[y] * Delta) + (r * value_matrix[y,-x])
        value_matrix[y,-x-1] = value_matrix[y,-x] - Theta * dt
        value_matrix[0,-x-1] = 2 * value_matrix[1,-x-1] - value_matrix[2,-x-1]
value_df = pd.DataFrame(value_matrix)
value_df = value_df.set_index(asset_price)
print(value_df[0])

plot_df = value_df.sort_index(ascending=True)
plot_df[0].plot()
plot_df[NTS-1].plot()
plt.show()
