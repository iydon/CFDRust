root = __import__('pathlib').Path(__file__).parent.absolute()
__import__('sys').path.insert(0, str(root.parents[1]))


import itertools as it
import re
import typing as t

import matplotlib.pyplot as plt
import numpy as np

from test import Fortran


def plot_solution(
    numerical: np.ndarray, analytical: np.ndarray, path: str,
    alpha: float = 0.5,
) -> None:
    x = np.linspace(-1.0, 1.0, numerical.size)
    fig, (ax1, ax2) = plt.subplots(
        nrows=2, ncols=1, sharex=True, figsize=(12, 6),
        gridspec_kw={'height_ratios': [3, 1]},
    )
    # numerical and analytical solutions
    ax1.plot(x, numerical, alpha=alpha, label='Numerical Solution')
    ax1.plot(x, analytical, alpha=alpha, label='Analytical Solution')
    ax1.grid()
    ax1.legend()
    # error between solutions
    ax2.plot(x, numerical-analytical, label='Error (N-A)')
    ax2.grid()
    ax2.legend()
    fig.savefig(path, bbox_inches='tight', transparent=True)

def plot_dimensionless_time(Ns: t.List[int], nutLs: t.List[float], alpha: float = 0.5) -> None:
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(Ns, nutLs, 'o-', alpha=alpha)
    ax.grid()
    ax.set_xlabel('Grid Resolution (N)')
    ax.set_ylabel('Dimensionless Time (νt/L²)')
    fig.savefig('dimensionless_time.png', bbox_inches='tight', transparent=True)

def plot_error(norms: np.ndarray, Ns: np.ndarray, nutLs: np.ndarray) -> None:
    fig, ax = plt.subplots(1, figsize=(16, 9))
    im = ax.imshow(np.log10(norms))
    ax.set_title(r'$\log_{10}||\mathrm{numerical}-\mathrm{analytical}||_2$')
    ax.set_xticks(np.arange(Ns.size), rotation=90, labels=Ns)
    ax.set_yticks(np.arange(nutLs.size), labels=map('{:.6f}'.format, nutLs))
    ax.set_xlabel('Grid Resolution (N)')
    ax.set_ylabel('Dimensionless Time (νt/L²)')
    fig.colorbar(im, ax=ax)
    fig.savefig('error.png', bbox_inches='tight', transparent=True)


if __name__ == '__main__':
    fortran = Fortran(root/'unidirectional_flow.f90').compile()
    # 1. plot and compare the velocity profiles from different grid resolutions
    for N, nutL in it.product([8, 16, 32], [0.2, 1.0, 10.0]):
        fortran.communicate(f'{N}\n{nutL}')
        numerical, analytical = np.loadtxt('numerical.dat'), np.loadtxt('analytical.dat')
        plot_solution(numerical, analytical, f'{N}-{nutL}.png')
    # 2. find out at what dimensionless time, νt/L², the velocity at the center of the channel (y=0) reaches the value of 0.99u₀
    Ns = list(range(2, 101, 2))
    nutLs = [None] * len(Ns)
    pattern = re.compile(r'\[N\][^\n]+\n')
    for ith, N in enumerate(Ns):
        stdout, _ = fortran.communicate(f'{N}\n0.0')
        nutLs[ith] = float(pattern.findall(stdout.decode())[0][33:43]) # hard coding
    plot_dimensionless_time(Ns, nutLs)
    # 3. use the analytical solution as the benchmark, plot the errors
    Ns = np.arange(8, 129, 2)
    nutLs = np.linspace(0.0, 2.0, Ns.size+1)[1:]
    norms = np.ndarray(shape=(nutLs.size, Ns.size))
    for (ith, nutL), (jth, N) in it.product(enumerate(nutLs), enumerate(Ns)):
        fortran.communicate(f'{N}\n{nutL}')
        numerical, analytical = np.loadtxt('numerical.dat'), np.loadtxt('analytical.dat')
        norms[ith, jth] = np.linalg.norm(numerical-analytical, ord=2)
    plot_error(norms, Ns, nutLs)
