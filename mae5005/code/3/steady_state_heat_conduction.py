root = __import__('pathlib').Path(__file__).parent.absolute()
__import__('sys').path.insert(0, str(root.parents[1]))


from cProfile import label
import itertools as it
import re
import typing as t

import matplotlib.pyplot as plt
import numpy as np

from test import Fortran


def plot_contour(numerical: np.ndarray, path: str) -> None:
    ny, nx = numerical.shape
    fig, ax = plt.subplots(1, figsize=(16, 9))
    im = ax.imshow(numerical)
    co = ax.contour(numerical, levels=10, colors='black')
    ax.clabel(co, inline=True, fontsize=8)
    ax.set_title('$T/T_0$')
    ax.set_xlabel(f'${nx}*x/L$')
    ax.set_ylabel(f'${ny}*y/L$')
    fig.colorbar(im, ax=ax)
    fig.savefig(path, bbox_inches='tight', transparent=True)

def plot_center(numerical: np.ndarray, path: str) -> None:
    ny, nx = numerical.shape
    assert nx == ny
    fig, ax = plt.subplots(1, figsize=(16, 9))
    ax.plot(numerical[:, (nx-1)//2], label='x/L = 0.5')
    ax.plot(numerical[(ny-1)//2, :], label='y/L = 0.5')
    ax.set_title('Temperature profiles at the two center lines')
    ax.set_xlabel(f'${ny}*y/L$ or ${nx}*x/L$')
    ax.set_ylabel(f'$T/T_0$')
    ax.grid()
    ax.legend()
    fig.savefig(path, bbox_inches='tight', transparent=True)


if __name__ == '__main__':
    nx = ny = 99

    fortran = Fortran(root/'steady_state_heat_conduction.f90').compile()
    fortran.communicate(f'{nx}\n{ny}')
    numerical = np.loadtxt('numerical.dat').reshape((ny+1, nx+1))
    plot_contour(numerical, 'contour.png')
    plot_center(numerical, 'center.png')
