import matplotlib.pyplot as plt
import numpy as np

from numpy import pi


a = 1
C = 0.5
x = np.linspace(0, 2, 256)

for N in [20, 160]:
    fig, ax = plt.subplots(1, figsize=(16, 9))
    for t, color in [(1.0, '#5EAA5F'), (2.0, '#9D6AB9')]:
        dx = 2 / N
        dt = C * dx / a
        # analytical solution
        ua = np.sin(pi*(x-a*t)) + 0.2*np.sin(4*pi*(x-a*t))
        # modified analytical solution
        b2 = a*dx/2 * (1-C)
        b3 = -a*dx**2/6 * (1-3*C+2*C**2)
        b4 = a*dx**3/24 * (1-7*C+12*C**2-6*C**3)
        c1 = -b2*pi**2 + b4*pi**4
        c2 = -a*pi - b3*pi**3
        c3 = -b2*(4*pi)**2 + b4*(4*pi)**4
        c4 = -4*a*pi - b3*(4*pi)**3
        um = np.exp(c1*t)*np.sin(pi*x+c2*t) + 0.2*np.exp(c3*t)*np.sin(4*pi*x+c4*t)
        # numerical solution
        xn = (np.arange(N+1) - 0.5) * dx
        un = np.sin(pi*xn) + 0.2*np.sin(4*pi*xn)
        for _ in range(int(t/dt)):
            un[1:] = un[1:] - C*(un[1:]-un[:-1])
            un[0] = un[-1]
        xn, un = xn[1:], un[1:]
        # image
        ax.plot(x, ua, color=color, linestyle='-', label=f't={t}, $u$')
        ax.plot(x, um, color=color, linestyle='--', label=f't={t}, $u_M$')
        ax.plot(xn, un, color=color, marker='o', linewidth=0, label=f't={t}, $u_N$')
    ax.set_title(f'N = {N}  (Upwind Scheme)')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u(x, t)$')
    ax.grid()
    ax.legend()
    fig.savefig(f'upwind-{N}.png', bbox_inches='tight', transparent=True)
