import matplotlib.pyplot as plt
import numpy as np

from numpy import pi


def rmse(u_true, u_pred):
    return np.sqrt(np.mean((u_true - u_pred) ** 2))


a = 1
C = 0.5
t = 1.0
Ns = np.arange(20, 160+2, 2)

schemes = {'Upwind': Ns.astype(float), 'Lax-Wendroff': Ns.astype(float)}
for ith, N in enumerate(Ns):
    dx = 2 / N
    dt = C * dx / a
    x = np.linspace(0, 2, N)
    # analytical
    ua = np.sin(pi*(x-a*t)) + 0.2*np.sin(4*pi*(x-a*t))
    # upwind
    b2 = a*dx/2 * (1-C)
    b3 = -a*dx**2/6 * (1-3*C+2*C**2)
    b4 = a*dx**3/24 * (1-7*C+12*C**2-6*C**3)
    c1 = -b2*pi**2 + b4*pi**4
    c2 = -a*pi - b3*pi**3
    c3 = -b2*(4*pi)**2 + b4*(4*pi)**4
    c4 = -4*a*pi - b3*(4*pi)**3
    uu = np.exp(c1*t)*np.sin(pi*x+c2*t) + 0.2*np.exp(c3*t)*np.sin(4*pi*x+c4*t)
    # lax_wendroff
    b3 = -a*dx**2/6 * (1-C**2)
    b4 = -a*dx**3/8 * C * (1-C**2)
    c1 = b4*pi**4
    c2 = -a*pi - b3*pi**3
    c3 = b4*(4*pi)**4
    c4 = -4*a*pi - b3*(4*pi)**3
    ul = np.exp(c1*t)*np.sin(pi*x+c2*t) + 0.2*np.exp(c3*t)*np.sin(4*pi*x+c4*t)
    # errors
    schemes['Upwind'][ith] = rmse(ua, uu)
    schemes['Lax-Wendroff'][ith] = rmse(ua, ul)
# image
fig, ax = plt.subplots(1, figsize=(16, 9))
for scheme, errors in schemes.items():
    ax.plot(Ns, errors, label=scheme)
    for N in [20, 40, 80, 160]:
        error = errors[Ns == N][0]
        ax.scatter([N], [error], c='gray')
        ax.text(N, error, f'{error:.7}', ha='center')
ax.set_title('L2 Errors Predicted Using the Modified Equation')
ax.set_xlabel('$N$')
ax.set_ylabel('L2 Errors')
ax.set_ylim([0.0, 0.23])
ax.grid()
ax.legend()
fig.savefig(f'modified.png', bbox_inches='tight', transparent=True)
