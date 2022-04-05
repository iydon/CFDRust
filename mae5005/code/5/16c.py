import matplotlib.pyplot as plt
import numpy as np

from numpy import pi


def rmse(u_true, u_pred):
    return np.sqrt(np.mean((u_true - u_pred) ** 2))


a = 1
T = 5.0
N = 80
dx = 2 / N
Cs = [0.5, 1.0, 1.2]

fig, ax = plt.subplots(1, figsize=(16, 9))
xa = np.linspace(0, 2, N)
for C in Cs:
    dt = C * dx / a
    nt = int(T / dt)
    times = np.linspace(0, T, nt)
    # upwind
    xu = (np.arange(N+1) - 0.5) * dx
    uu = np.sin(pi*xu) + 0.2*np.sin(4*pi*xu)
    # lax_wendrof
    xl = (np.arange(N+2) - 0.5) * dx
    ul = np.sin(pi*xl) + 0.2*np.sin(4*pi*xl)
    # errors
    error_u = np.zeros(nt)
    error_l = np.zeros(nt)
    for ith in range(nt):
        t = (ith+1) * dt
        # analytical
        ua = np.sin(pi*(xa-a*t)) + 0.2*np.sin(4*pi*(xa-a*t))
        # upwind
        uu[1:] = uu[1:] - C*(uu[1:]-uu[:-1])
        uu[0] = uu[-1]
        # lax_wendroff
        ul[1:-1] = ul[1:-1] - C/2*(ul[2:]-ul[:-2]) + C**2/2*(ul[2:]-2*ul[1:-1]+ul[:-2])
        ul[0], ul[-1] = ul[-2], ul[1]
        # errors
        error_u[ith] = rmse(ua, uu[1:])
        error_l[ith] = rmse(ua, ul[1:-1])
    ax.semilogy(times, error_u, linestyle='-', label=f'C = {C} (Upwind Scheme)')
    ax.semilogy(times, error_l, linestyle='--', label=f'C = {C} (Lax-Wendroff Scheme)')
ax.set_xlabel('Time')
ax.set_ylabel('L2 Errors')
ax.set_ylim([1e-2, 1e2])
ax.grid()
ax.legend()
fig.savefig(f'time.png', bbox_inches='tight', transparent=True)
