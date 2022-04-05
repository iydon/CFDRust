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
    # analytical
    xa = np.linspace(0, 2, N)
    ua = np.sin(pi*(xa-a*t)) + 0.2*np.sin(4*pi*(xa-a*t))
    # upwind
    xu = (np.arange(N+1) - 0.5) * dx
    uu = np.sin(pi*xu) + 0.2*np.sin(4*pi*xu)
    for _ in range(int(t/dt)):
        uu[1:] = uu[1:] - C*(uu[1:]-uu[:-1])
        uu[0] = uu[-1]
    xu, uu = xu[1:], uu[1:]
    # lax_wendroff
    xl = (np.arange(N+2) - 0.5) * dx
    ul = np.sin(pi*xl) + 0.2*np.sin(4*pi*xl)
    for _ in range(int(t/dt)):
        ul[1:-1] = ul[1:-1] - C/2*(ul[2:]-ul[:-2]) + C**2/2*(ul[2:]-2*ul[1:-1]+ul[:-2])
        ul[0], ul[-1] = ul[-2], ul[1]
    xl, ul = xl[1:-1], ul[1:-1]
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
ax.set_title('L2 Errors')
ax.set_xlabel('$N$')
ax.set_ylabel('L2 Errors')
ax.set_ylim([0.0, 0.23])
ax.grid()
ax.legend()
fig.savefig(f'N.png', bbox_inches='tight', transparent=True)
