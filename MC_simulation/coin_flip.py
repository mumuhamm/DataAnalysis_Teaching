#Monte Carlon hands on - coin toss / cross check with binomial distribution
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

# 0 - heads 
# 1 - tails 
def trial():
    return random.randint(0,1)
print(trial())

list = []
binom_data = []
binom_prob = []
def experiment(n):
    results = 0
    for i in range(n):
        toss = trial()
        results = results + toss
        probability = results/(i+1)
        list.append(probability)
    return results/n

expt = experiment(1000)
print("Probability of either getting a heads or tails in 1000 experiments from MC simulation")
print(list)
#print("Probability of getting heads or tails :", expt)

k = 100 # trial
p = 0.5 # probability
size=1000
binom_data.append((np.random.binomial(k, p, size=size))/100)
print("Length of binomial data :", len(binom_data))
for i in range(0, len(binom_data)):
    binom_prob.append(binom_data[i])
    print(binom_data[i])


#=========================Plot
plt.axhline(y=0.50, color = 'r', linestyle = '-')
plt.xlabel("Trials")
plt.ylabel("Probability")
plt.plot(list,'o')
plt.plot(binom_prob,'o')
plt.show()
