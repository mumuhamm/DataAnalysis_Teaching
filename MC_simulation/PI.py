#import tk
#import turtle
#from turtle import *
# Turtle is not working , tkinter issue , possibly my os version has the origin
# value of pi , can be cross cheked with math library

import math
import random 
import matplotlib.pyplot as plt

count_pt_in_circle = 0
count_pt_out_circle = 0

pi_values = []

for i in range(5):
    for j in range(10000):
        x = random.randrange(-100,100)
        y = random.randrange(-100,100)
        if(x**2+y**2 > 100**2):
          count_pt_out_circle = count_pt_out_circle + 1
        else:
            count_pt_in_circle = count_pt_in_circle + 1
        pi = 4.0 * count_pt_in_circle /( count_pt_out_circle + count_pt_in_circle)
        pi_values.append(pi)
        error_pi = [abs(math.pi-pi) for pi in pi_values]
        
    print(pi_values[-1])

plot1 = plt.figure(1)
plt.axhline(y=math.pi, color = 'g', linestyle = 'dashed')
plt.plot(pi_values)
plt.xlabel("Iteration")
plt.ylabel("$\pi$ value")
plot2 = plt.figure(2)
plt.axhline(y=0, color = 'r', linestyle = 'dashdot')
plt.plot(error_pi)
plt.xlabel("Iteration")
plt.ylabel("$\Delta\pi$ : Error")
plt.show()
