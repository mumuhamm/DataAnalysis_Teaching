import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

"""
This is all from the tutorial located at :
http://scipy-lectures.github.io/intro/matplotlib/matplotlib.html
"""

pl.figure(figsize=(10, 6), dpi=80)
pl.subplot(1, 1, 1)
X = np.linspace(-5, 5, 500, endpoint=True)
C = (1/X**2)-5
P = X - X - 0.1

pl.xlim(X.min() * 1.1, X.max() * 1.1)
pl.ylim(C.min() * 1.1, C.max() * 1.1)

"""
Alters the position of the axis - moves them to the centre
"""
ax = pl.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

pl.plot(X, C, color="blue",  linewidth=4, linestyle="-", 
              label="y = 4 - 1/x^2")

pl.legend(loc='upper left')


for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
    label.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.65))


ylim = ax.get_ylim()
plt.ylim((-7,20))
plt.vlines(3, ylim[0], ylim[1])
plt.show()
