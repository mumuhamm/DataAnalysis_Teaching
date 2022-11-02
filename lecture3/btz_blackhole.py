import numpy as np
from mpmath import hyp2f1, mpc, mpf, gamma, mp, sin, cos, sqrt, exp, log, pi, conj
import plotly.graph_objects as go
import matplotlib.pyplot as plt

class BTZ:
    # class object determining BTZ black holes

    def __init__(self, M, r_plus, mu2=-0.75):
        self.M = M            # mass of the black hole
        self.r_plus = r_plus  # square root of the mass of the black hole
        self.mu2 = mu2        # effective mass for the massless conf. coup. scalar field

    def g00(self, z):
        # coefficient of dt**2 of the metric tensor
        return self.M*z/(z-1)

    def omega_tilde_00(self, omega, z):
        # blueshift effect on omega
        return mp.sqrt(abs(self.g00(z)))*abs(omega)

    def a(self, l, omega):
        return 1/2 * (1 + sqrt(1 + self.mu2) + 1j*(omega + l)/self.r_plus)
    # parameter a of the hypergeometric solutions

    def b(self, l, omega):
        return 1/2 * (1 + sqrt(1 + self.mu2) + 1j*(omega - l)/self.r_plus)
    # parameter b of the hypergeometric solutions

    def c(self, omega):
        return 1 + 1j*omega/self.r_plus
    # parameter c of the hypergeometric solutions

    def A(self, l, omega):
        return mp.gamma(self.c(omega)) * \
    mp.gamma(self.c(omega) - self.a(l, omega) - self.b(l, omega)) / \
    mp.gamma(self.c(omega) - self.a(l, omega))/mp.gamma(self.c(omega) - self.b(l, omega))
    # parameter A of the normalization constant of the transition rate

    def B(self, l, omega):
        # parameter B of the normalization constant of the transition rate
        return mp.gamma(self.c(omega))*mp.gamma(self.a(l, omega) +\
                                                    self.b(l, omega) - self.c(omega)) /\
            gamma(self.a(l, omega))/mp.gamma(self.b(l, omega))

    def factor_C(self):
        # parameter C of the normalization constant of the transition rate
        return 2/(1j*mp.pi*mp.sqrt(1 + self.mu2))

    def normalization(self, gamma, l, omega):
        return self.factor_C()*(conj(self.A(l, omega))*self.B(l, omega) -
                                self.A(l, omega) * conj(self.B(l, omega))) /\
            (abs(self.B(l, omega)*cos(gamma) - self.A(l, omega)*sin(gamma))**2)
    # normalization parameter that goes in the transition rate

    def alpha(self, omega):
        return 1j * omega / (2 * self.r_plus)
    # parameter of the ansatz of the radial solutions

    def beta_plus(self):
        # parameter of the ansatz of the radial solutions,  equals beta
        return 1/2 + mp.sqrt(1 + self.mu2)/2

    def beta_minus(self):
        # parameter of the ansatz of the radial solutions,  equals 1 - beta
        return 1/2 - mp.sqrt(1 + self.mu2)/2

    def R11(self, l, gamma, omega, z):
        # R11 and R22 form a basis of radial solutions at AdS boundary
        return z**self.alpha(omega)*(1-z)**self.beta_plus() * \
            hyp2f1(self.a(l, omega), self.b(l, omega), self.a(l, omega) +
                   self.b(l, omega) - self.c(omega) + 1, 1 - z)

    def R21(self, l, gamma, omega, z):
        # R11 and R22 form a basis of radial solutions at AdS boundary
        return z**self.alpha(omega)*(1-z)**self.beta_minus() * \
            hyp2f1(self.c(omega) - self.a(l, omega), self.c(omega) - self.b(l, omega),
                   self.c(omega) - self.a(l, omega) - self.b(l, omega) + 1, 1 - z)

    def R(self, l, gamma, omega, z):
        return mp.cos(gamma)*self.R11(l, gamma, omega, z) + mp.sin(gamma) *\
            self.R21(l, gamma, omega, z)
    # radial solution at the AdS boundary compatible with Robin b.c.

    def z_of_TH(self, TH):
        return 1/(1 + 4*pi**2*TH**2)
    # detector's position as function of the Hawking temperature

    def TH_of_z(self, z):
        return self.r_plus/(2*pi*sqrt(abs(self.g00(z))))
    # Hawking temperature as function of detector's position

    def transition_rate_0(self, gamma, l, Egap, z):
        # l-term constribution to the transition rate
        # for the ground state on the static BTZ black hole
        # as a function of the boundary contidion parametrized by gamma,
        # the quantum number l, the energy gap Egap and the position
        # of the detector z.
        if Egap < 0:
            x = self.normalization(gamma, l, self.omega_tilde_00(Egap, z)) * \
                self.R(l, gamma, self.omega_tilde_00(Egap, z), z)**2
            return x.real, x.imag
        else:
            return 0, 0
    def transition_rate_T(self, gamma, l, Egap, z):
        # l-term constribution to the transition rate
        # for the KMS state on the static BTZ black hole
        # as a function of the boundary contidion parametrized by gamma,
        # the quantum number l, the energy gap Egap and the position
        # of the detector z.
        if Egap < 0:
            x = 1/(1 - exp(-self.omega_tilde_00(Egap, z)/self.TH_of_z(z))) * \
                self.normalization(gamma, l, self.omega_tilde_00(Egap, z)) * \
                self.R(l, gamma, self.omega_tilde_00(Egap, z), z)**2
            return x.real, x.imag
        else:
            x = 1/(exp(self.omega_tilde_00(Egap, z)/self.TH_of_z(z)) - 1) * \
                self.normalization(gamma, l, self.omega_tilde_00(Egap, z)) * \
                self.R(l, gamma, self.omega_tilde_00(Egap, z), z)**2
            return x.real, x.imag
            
btz1 = BTZ(1, 1)
normalization_factor = 2.75

N = 200
x = np.linspace(1e-5, 1-1e-5, N)

lp, Egapp = 0, -0.1

gamma0 = 0
def func0(p):
    return float(btz1.transition_rate_0(gamma0,lp,Egapp,p)[0])/normalization_factor
vfunc0 = np.vectorize(func0)
y0 = (vfunc0((x)))

gamma1 = 0.25*np.pi
def func1(p):
    return float(btz1.transition_rate_0(gamma1,lp,Egapp,p)[0])/normalization_factor
vfunc1 = np.vectorize(func1)
y1 = (vfunc1((x)))

gamma2 = 0.40*np.pi
def func2(p):
    return float(btz1.transition_rate_0(gamma2,lp,Egapp,p)[0])/normalization_factor
vfunc2 = np.vectorize(func2)
y2 = (vfunc2((x)))

gamma3 = 0.47*np.pi
def func3(p):
    return float(btz1.transition_rate_0(gamma3,lp,Egapp,p)[0])/normalization_factor
vfunc3 = np.vectorize(func3)
y3 = (vfunc3((x)))

gamma4 = 0.5*np.pi
def func4(p):
    return float(btz1.transition_rate_0(gamma4,lp,Egapp,p)[0])/normalization_factor
vfunc4 = np.vectorize(func4)
y4 = (vfunc4((x)))





fig, ax = plt.subplots()
ax.plot(x, y4, '-b', label='Newmann')
ax.plot(x, y3, '--r', label='Robin 3')
ax.plot(x, y2, '-.k', label='Robin 2')
ax.plot(x, y1, ':g', label='Robin 1')
ax.plot(x, y0, '-.y', label='Dirichlet')
ax.axis('equal')
leg = ax.legend()
ax.legend(loc='upper left', frameon=False)
plt.xlabel('$z$')
plt.ylabel('$\dot{\mathcal{F}}$')
plt.show()

'''plt.scatter(x=x,y=y4,color="k", name='Newmann' )
plt.scatter(x=x,y=y3,color="r", name='Newmann')
plt.scatter(x=x,y=y2,color="b", name='Newmann')
plt.scatter(x=x,y=y1,color="y", name='Newmann')
plt.scatter(x=x,y=y0,color="g", name='Newmann')



fig.update_layout(
    title="Transition rate for the Ground state as a function of position z on the static BTZ black hole",
    xaxis_title=r"$z$",
    yaxis_title=r"$\dot{\mathcal{F}}$",
    legend_title="Boundary Conditions"
#     ,font=dict(
#         family="Helvetica",
#         size=18,
#         color="RebeccaPurple"
#     )
)

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y4, name='Newmann'))
fig.add_trace(go.Scatter(x=x, y=y3, name='Robin 3'))
fig.add_trace(go.Scatter(x=x, y=y2, name='Robin 2'))
fig.add_trace(go.Scatter(x=x, y=y1, name='Robin 1'))
fig.add_trace(go.Scatter(x=x, y=y0, name='Dirichlet'))


fig.update_layout(
    title="Transition rate for the Ground state as a function of position z on the static BTZ black hole",
    xaxis_title=r"$z$",
    yaxis_title=r"$\dot{\mathcal{F}}$",
    legend_title="Boundary Conditions"
#     ,font=dict(
#         family="Helvetica",
#         size=18,
#         color="RebeccaPurple"
#     )
)

fig.show()'''
