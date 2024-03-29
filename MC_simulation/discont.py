import matplotlib.pyplot as plt 
import numpy as np
from sympy import sympify, lambdify
from sympy.abc import x

fig = plt.figure(1) 
ax = fig.add_subplot(111) 

# set up axis 
ax.spines['left'].set_position('zero') 
ax.spines['right'].set_color('none') 
ax.spines['bottom'].set_position('zero') 
ax.spines['top'].set_color('none') 
ax.xaxis.set_ticks_position('bottom') 
ax.yaxis.set_ticks_position('left') 

# setup x and y ranges and precision
xx = np.arange(-0.5,5.5,0.01) 

# draw my curve 
myfunction=sympify(1/(x-2))
mylambdifiedfunction=lambdify(x,myfunction,'numpy')
ax.plot(xx, mylambdifiedfunction(xx),zorder=100,linewidth=3,color='red') 

#set bounds 
ax.set_xbound(-1,6)
ax.set_ybound(-4,4) 

plt.show()
