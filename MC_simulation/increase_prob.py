import random
import matplotlib.pyplot as plt
import math 
import numpy as np

#four choice, the player is allowed to switch the choice after the first random choice, the idea is here to look
#whether the next choice (switching the choice) will increase thepossibility of win or not
#x1 is the desired choice  , prob 1/4
#x2 is not the desired choice, prob 3/4
# three x2 and one x1
# Please extend it to n dimension 

choice = ["x1", "x2", "x2", "x2"]
switch_prob = []
not_switch_prob = []

def experiment(n):
    switch = 0
    not_switch = 0
    for i in range(n):
        random.shuffle(choice)
        j = random.randrange(3)
        if choice[j] !='x1':
           switch = switch +1
        else:
            not_switch = not_switch + 1
        switch_prob.append(switch/(i+1))
        not_switch_prob.append(not_switch/(i+1))
        
    print('Probability with switch choice : ', switch_prob[-1])
    print('Probability with not switch choice : ', not_switch_prob[-1])
   
experiment(1000)
plt.axhline(y=0.75, color = 'g', linestyle = 'dashdot')
plt.axhline(y=0.25, color = 'g', linestyle = 'dashdot')
plt.plot(switch_prob, '--b', label='switch')
plt.plot(not_switch_prob, c='r', label='not switch')
plt.xlabel('Iteration : $n$')
plt.ylabel('Probability : $\mathcal{P}(x)$')
plt.legend(loc='upper right')
plt.show()
