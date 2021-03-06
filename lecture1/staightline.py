import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

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
  
