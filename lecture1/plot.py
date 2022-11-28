import matplotlib.pyplot as plt 

def plotx(axs, bxs, cys, dys, dyserr):
    plt.xlabel('x', fontsize=15)
    plt.ylabel('f(x)', fontsize=15)
    plt.plot(axs, cys, 'r-', label='quadratic')
    plt.errorbar(bxs, dys, dyserr, fmt='b:D', label='some function with error')
    plt.legend()
    plt.show()

axs = [0.1*i for i in range (60)]
bxs = [i for i in range (7)]
cys = [x**2 for x in axs]
dys = [x+x**2 for x in bxs]
dyserr = [0.1*y for y in dys]

plotx(axs, bxs, cys, dys, dyserr)
