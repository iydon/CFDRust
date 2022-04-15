import collections as c

import matplotlib.pyplot as plt
import numpy as np


def minmod(a, *bs, func=np.min):
    ans = np.zeros_like(a)
    sign = np.sign(a)
    flag = np.all([sign==np.sign(b) for b in bs], axis=0)
    ans[flag] = sign[flag] * func([np.abs(b[flag]) for b in (a,)+bs], axis=0)
    return ans

def maxmod(a, *bs):
    return minmod(a, *bs, func=np.max)

def total_variation(u):
    return np.sum(np.abs(np.diff(u)))


a, dt, dx = 1.0, 0.05, 0.1
cfl = a * dt / dx
nts = {1, 4, 20, 40, 75}

x = (np.arange(-20.0, 40.0) + 0.5) * dx
u = lambda t: 2.0 * (x<a*t)
u0 = u(0.0)
data = c.defaultdict(dict)
tvs = c.defaultdict(dict)

# downwind slope
ud = np.zeros(u0.size+2)
ud[1:-1] = u0
for nt in range(1, max(nts)+1):
    ud[[0, -1]] = ud[[1, -2]]
    sigma1 = (ud[2:] - ud[1:-1]) / dx
    sigma2 = (ud[1:-1] - ud[:-2]) / dx
    ud[1:-1] = ud[1:-1] - cfl*(ud[1:-1]-ud[:-2]) - cfl/2.0*(dx-a*dt)*(sigma1-sigma2)
    if nt in nts:
        data[nt]['downwind slope'] = ud[1:-1].copy()
    tvs['downwind slope'][nt] = total_variation(ud[1:-1])
# van Leer limiter
uv = np.zeros(u0.size+4)
uv[2:-2] = u0
for nt in range(1, max(nts)+1):
    uv[[0, 1, -2, -1]] = uv[[2, 2, -3, -3]]
    sigma1 = minmod(2.0*(uv[3:-1]-uv[2:-2])/dx, (uv[3:-1]-uv[1:-3])/2.0/dx, 2.0*(uv[2:-2]-uv[1:-3])/dx)
    sigma2 = minmod(2.0*(uv[2:-2]-uv[1:-3])/dx, (uv[2:-2]-uv[:-4])/2.0/dx, 2.0*(uv[2:-2]-uv[:-4])/dx)
    uv[2:-2] = uv[2:-2] - cfl*(uv[2:-2]-uv[1:-3]) - cfl/2.0*(dx-a*dt)*(sigma1-sigma2)
    if nt in nts:
        data[nt]['van-Leer slope limiter'] = uv[2:-2].copy()
    tvs['van-Leer slope limiter'][nt] = total_variation(uv[2:-2])
# SUPERBEE limiter
us = np.zeros(u0.size+4)
us[2:-2] = u0
for nt in range(1, max(nts)+1):
    us[[0, 1, -2, -1]] = us[[2, 2, -3, -3]]
    sigma1 = maxmod(
        minmod(2.0*(us[2:-2]-us[1:-3])/dx, (us[3:-1]-us[2:-2])/dx),
        minmod((us[2:-2]-us[1:-3])/dx, 2.0*(us[3:-1]-us[2:-2])/dx),
    )
    sigma2 = maxmod(
        minmod(2.0*(us[1:-3]-us[:-4])/dx, (us[2:-2]-us[1:-3])/dx),
        minmod((us[1:-3]-us[:-4])/dx, 2.0*(us[2:-2]-us[1:-3])/dx),
    )
    us[2:-2] = us[2:-2] - cfl*(us[2:-2]-us[1:-3]) - cfl/2.0*(dx-a*dt)*(sigma1-sigma2)
    if nt in nts:
        data[nt]['SUPERBEE slope limiter'] = us[2:-2].copy()
    tvs['SUPERBEE slope limiter'][nt] = total_variation(us[2:-2])

# figure
for nt, values in data.items():
    fig, ax = plt.subplots(1, figsize=(16, 9))
    ax.plot(x, u(nt*dt), label='analytical solution')
    for label, y in values.items():
        ax.plot(x, y, label=label)
    ax.set_ylim([-0.1, 2.4])
    ax.set_title(f'1D unsteady advection problem (timestep={nt})')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u$')
    ax.grid()
    ax.legend()
    fig.savefig(f'1d_unsteady_advection_problem-{nt}.png', bbox_inches='tight', transparent=True)

for nt in range(1, max(nts)+1):
    tvs['analytical solution'][nt] = total_variation(u(nt*dt))
fig, ax = plt.subplots(1, figsize=(16, 9))
for label, values in tvs.items():
    ax.plot(values.keys(), values.values(), linestyle='--', label=label)
ax.set_xlabel('Step')
ax.set_ylabel('Total Variation (TV)')
ax.grid()
ax.legend()
fig.savefig('total_variation.png', bbox_inches='tight', transparent=True)
