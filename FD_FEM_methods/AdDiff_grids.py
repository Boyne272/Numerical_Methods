#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Just some code to animate the time evolution of the
advection diffusion problem with dirichlet BCs

Created on Mon Dec  3 21:34:09 2018
@author: Richard Boyne rmb115@ic.ac.uk
"""

import numpy as np
from my_animate import animate_array, percent_print
from Part1 import exact_solution
import matplotlib.pyplot as plt


def ghost_dis_mat(x, lims, kappa, U):
    "lims are the edges with bcs that will be ghost nodes"
    # find dx array
    tmp_x = np.array(lims[:1] + x.tolist() + lims[1:])
    dx = np.diff(tmp_x)
    N = len(x)

    # set up matrix
    A = np.zeros([N, N+2])
    tmp1, tmp2 = -2 * kappa, U / 2
    for i in range(1, N+1):
        sigma = (dx[i] + dx[i-1]) * dx[i] * dx[i-1]
        A[i-1, i-1] = -dx[i] * tmp1 / sigma - tmp2 * (-1/dx[i-1])
        A[i-1, i] = (dx[i] + dx[i-1]) * tmp1 / sigma - \
            tmp2 * ((1/dx[i-1]) - (1/dx[i]))
        A[i-1, i+1] = -dx[i-1] * tmp1 / sigma - tmp2 * (1/dx[i])

    # add in the boundary conditions
    a, b = np.zeros([N+2]), np.zeros([N+2])
    a[0] = 1.
    b[-1] = 1.
    A = np.vstack((a, A, b))

    return A


# %%
options = [0, 2, 'ani']
marker = 'o-'

# set up the system
if 0 in options:
    # parameters
    kappa, U = 0.01, .3
    t_st, t_max, dt = 0., 8, 1e-3
    lims = [0., 1.]
    x = np.arange(0.005, 1.0, 0.005)
    x = np.array(x)
    x_full = np.array(lims[:1] + x.tolist() + lims[1:])

    # constants
    times = np.arange(t_st, t_max + dt, dt)
    t_num = len(times)
    N = len(x)
    dx = np.diff(x).mean()  # excludes ghost nodes
    Ident = np.eye(N+2)
    A = ghost_dis_mat(x, lims, kappa, U)

    # stability
    cell_Pe = dx * U / kappa
    CFL = U * dt / dx
    r = CFL / cell_Pe
    print('cell_Pe < 2', cell_Pe < 2,
          '\nCFL < 1    ', CFL < 1,
          '\n|r| < 0.5  ', 0 < r < 0.5)
    print('Pe No. = ', lims[-1] * U / kappa)


# solve FTCS
if 2 in options:
    # set up the grid and initial conditions
    C = np.empty([t_num, N+2])
#    C[0, 1:N+1] = exact_solution(0.1, x, kappa, U)
    C[0, 0] = 1.
    C[0, N+1] = 0.

    # time step
    tmp = (dt * A + Ident)  # remove repeat operation
    tmp[0, 0], tmp[-1, -1] = 1., 1.  # ########### QUICK FIX!
    for i in range(t_num-1):
        C[i+1, :] = tmp @ C[i, :]
        percent_print(i, t_num)


# animate the results
if 'ani' in options:
    ani = animate_array(C[::20], x_full, dt)
    ani.times = times[::20]
    ani.interval = 20
#    ani.save = "wiggle_wiggle_wiggle.mp4"
    ani.animate()


# plot the results
if 'plt' in options:
    ax = plt.gca()
    n = int(len(times) / 5)
    for i, res in enumerate(C[::n]):
        ax.plot(x_full[::5], res[::5], marker, label='t=%.4f' % times[n*i])
    ax.legend()
    ax.set(xlabel='x', ylabel='C', title='Time Plots')
