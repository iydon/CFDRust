root = __import__('pathlib').Path(__file__).parent.absolute()
__import__('sys').path.insert(0, str(root.parents[1]))


import typing as t

import matplotlib.pyplot as plt
import numpy as np

from test import Fortran


def calc(alpha: float, nutL: float, N: int, L: float, nu: float) -> t.Tuple[int, float]:
    dt = 4.0 * alpha * L**2 / N**2 / nu
    nt = int(nutL * L**2 / nu / dt)
    return nt, dt


if __name__ == '__main__':
    fortran = Fortran(root/'1d_transient_laminar_flow.f90').compile()

    # (a)
    alpha, nutL = 0.32, 2.0
    N, L, nu = 16, 1.0, 0.1
    nt, dt = calc(alpha, nutL, N, L, nu)
    stdout, _ = fortran.communicate(f'{N}\n{nt}\n{L}\n{nu}\n{0.05}')
    print(repr(f'{N}\n{nt}\n{L}\n{nu}\n{dt}'))
    print(repr(stdout.decode()))
    exit()
    un, ua, um = map(np.loadtxt, ['numerical.dat', 'analytical.dat', 'modified.dat'])
    times = nu * (np.arange(nt)+1) * dt / L**2
    fig, ax = plt.subplots(1, figsize=(16, 9))
    ax.scatter(times, np.abs((um-ua)[:, N//2+1]), s=1, label=r'$|u_m-u|/u_0$')
    ax.scatter(times, np.abs((un-ua)[:, N//2+1]), s=1, label=r'$|u_N-u|/u_0$')
    ax.set_title('Temperature profiles at the two center lines')
    ax.set_xlabel('Dimensionless Time (νt/L²)')
    ax.set_ylabel('Normalized Error')
    ax.grid()
    ax.legend()
    fig.savefig('error.png', bbox_inches='tight', transparent=True)
    # (b)
    alpha = 0.32
    L, nu = 1.0, 0.1
    for nutL in [0.2, 10.0]:
        fig, ax = plt.subplots(1, figsize=(16, 9))
        for N in [8, 16, 32]:
            nt, dt = calc(alpha, nutL, N, L, nu)
            fortran.communicate(f'{N}\n{nt}\n{L}\n{nu}\n{dt}')
            un, ua, um = map(np.loadtxt, ['numerical.dat', 'analytical.dat', 'modified.dat'])
            ys = np.linspace(-L, L, N+1)
            ax.plot(ys, np.abs((um-ua)[-1, :]), '--', label=f'N={N} (analytical)')
            ax.scatter(ys, np.abs((un-ua)[-1, :]), label=f'N={N} (numerical)')
        ax.set_xlabel('y/L')
        ax.set_ylabel('Normalized Error')
        ax.set_title(f'α={alpha}, νt/L²={nutL}')
        ax.grid()
        ax.legend()
        fig.savefig(f'nutL-{nutL}.png', bbox_inches='tight', transparent=True)
    # (c)
    nutL = 2.0
    N, L, nu = 32, 1.0, 0.1
    fig, ax = plt.subplots(1, figsize=(16, 9))
    for alpha in [0.5, 0.505, 0.51, 0.52]:
        nt, dt = calc(alpha, nutL, N, L, nu)
        fortran.communicate(f'{N}\n{nt}\n{L}\n{nu}\n{dt}')
        un = np.loadtxt('numerical.dat')
        times = nu * (np.arange(nt)+1) * dt / L**2
        ax.scatter(times, un[:, N//2+1], s=1, label=f'α={alpha}')
    ax.set_ylim([0.0, 2.0])
    ax.set_xlabel('Dimensionless Time (νt/L²)')
    ax.set_ylabel('Velocity at Center-Line (y=0)')
    ax.grid()
    ax.legend()
    fig.savefig('alpha.png', bbox_inches='tight', transparent=True)
