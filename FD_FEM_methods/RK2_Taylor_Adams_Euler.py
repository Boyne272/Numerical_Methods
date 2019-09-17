#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collection of all code from tasks 1A, 1B, 1C, 1D for calling in the Report
Notebook

Created on Wed Nov 21 17:30:53 2018
@author: Richard Boyne rmb115@ic.ac.uk
"""

import numpy as np


# task-1A


def RK2_alpha(f, y0, t0, t_max, dt, alpha=0.5):
    """
    Integrates the given function by the RK2 method with the alpha value
    given.

    PARAMETERS
    ----------
    f : callable with inputs f(t, y)
        derivative function to be integrated. Must take the current time, t
        (float) and current solution values (1D array) and return the
        derivatives at those points in an array of the same size as
        the input array y.

    y0 : 1D array_like
        initial conditions, must be the same size as that parameters for
        function f. Must be indexable (if single variable problem place
        in list by itself)

    t0, t_max : float_like
        the range of times [t0, t_max] to integrate over with t_max > t0.
        Note that the final solution value returned will be at a
        time >= t_max depending on the choice of dt

    dt : float
        step size for the integration, for evlation at t_max must satisfy

            dt = (t_max - t0) / n       for integer n

    alpha : float_like, optional
        the alpha parameter for the RK2 equations, must be positive

    RETURNS
    -------

    times : 1D array
        an array of times where the solution was evaluated at

    y : nd array
        an array of shape (len(times), len(y0)) with the solutions at the
        time values given.
    """

    # validate inputs
    assert callable(f), 'f must be a function'
    assert alpha > 0, 'alpha must be positive'
    y_init = np.array(y0)
    assert y_init.shape != (), 'y0 must be iterable (list, array or tuple)'
    t0, t_max, dt, alpha = float(t0), float(t_max), float(dt), float(alpha)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0] = y_init

    # solve for every t
    for n, t in enumerate(times[:-1]):
        k1 = dt * f(t, y[n])
        k2 = dt * f(t + alpha*dt, y[n] + alpha*k1)
        y[n+1] = y[n] + k1 * (1 - 1/(2*alpha)) + k2 / (2*alpha)

    return times, y


def Euler_Forwards(f, y0, t0, t_max, dt):
    """
    Integrates the given function by the Forward Euler method.

    PARAMETERS
    ----------
    f : callable with inputs f(t, y)
        derivative function to be integrated. Must take the current time, t
        (float) and current solution values (1D array) and return the
        derivatives at those points in an array of the same size as
        the input array y.

    y0 : 1D array_like
        initial conditions, must be the same size as that parameters for
        function f. Must be indexable (if single variable problem place
        in list by itself)

    t0, t_max : float_like
        the range of times [t0, t_max] to integrate over with t_max > t0.
        Note that the final solution value returned will be at a
        time >= t_max depending on the choice of dt

    dt : float
        step size for the integration, for evlation at t_max must satisfy

            dt = (t_max - t0) / n       for integer n

    alpha : float_like, optional
        the alpha parameter for the RK2 equations, must be positive

    RETURNS
    -------

    times : 1D array
        an array of times where the solution was evaluated at

    y : nd array
        an array of shape (len(times), len(y0)) with the solutions at the
        time values given.
    """

    # validate inputs
    assert callable(f), 'f must be a function'
    y0 = np.array(y0)
    assert y0.shape != (), 'y0 must be iterable'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)
    y = np.zeros([times.size, y0.size])
    y[0] = y0

    for n, t in enumerate(times[:-1]):
        y[n+1] = y[n] + dt * f(t, y[n])

    return times, y


def Euler_Improved(f, y0, t0, t_max, dt):
    """
    Integrates the given function by the predicator corrector Euler method.

    PARAMETERS
    ----------
    f : callable with inputs f(t, y)
        derivative function to be integrated. Must take the current time, t
        (float) and current solution values (1D array) and return the
        derivatives at those points in an array of the same size as
        the input array y.

    y0 : 1D array_like
        initial conditions, must be the same size as that parameters for
        function f. Must be indexable (if single variable problem place
        in list by itself)

    t0, t_max : float_like
        the range of times [t0, t_max] to integrate over with t_max > t0.
        Note that the final solution value returned will be at a
        time >= t_max depending on the choice of dt

    dt : float
        step size for the integration, for evlation at t_max must satisfy

            dt = (t_max - t0) / n       for integer n

    alpha : float_like, optional
        the alpha parameter for the RK2 equations, must be positive

    RETURNS
    -------

    times : 1D array
        an array of times where the solution was evaluated at

    y : nd array
        an array of shape (len(times), len(y0)) with the solutions at the
        time values given.
    """

    # validate inputs
    assert callable(f), 'f must be a function'
    y0 = np.array(y0)
    assert y0.shape != (), 'y0 must be iterable'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)
    y = np.zeros([times.size, y0.size])
    y[0] = y0

    # solve
    for n, t in enumerate(times[:-1]):
        f_yn = f(t, y[n])
        y_tmp = y[n] + dt * f_yn
        y[n+1] = y[n] + 0.5 * dt * (f_yn + f(t+dt, y_tmp))

    return times, y


def correct_final(y, times, t_eval):
    """
    Correct the final value if wanted by linear interpolation between
    the final two points
    """
    dt = times[-1] - times[-2]
    grad = (y[-1] - y[-2]) / dt
    dx = t_eval - times[-2]
    y_final = y[-2] + dx * grad  # y value at t_max

    y[-1] = y_final
    times[-1] = t_eval
    return times, y


# Task 1B


def simple_taylor(f_0, f_1, f_2,
                  y0, t0, t_max, dt, exact_end=False):
    """
    Calculates the taylor series up to the first 2nd, 3rd and 4th terms using
    the total derivatives suppiled
    """

    # validate inputs
    assert all([callable(func) for func in (f_0, f_1, f_2)]), \
        'all functions must be callable'
    y_init = np.array(y0)
    assert y_init.shape != (), 'y0 must be iterable (list, array or tuple)'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y_0 = np.zeros([times.size, y_init.size])  # to hold solutions in
    y_1 = y_0.copy()
    y_2 = y_0.copy()
    y_0[0], y_1[0], y_2[0] = y_init, y_init, y_init

    # solve for every t
    for n, t in enumerate(times[:-1]):
        # solution to the 1st order
        y_0[n+1] = y_0[n] + f_0(t, y_0[n]) * dt
        # solution to the 2nd order
        y_1[n+1] = y_1[n] + f_0(t, y_1[n]) * dt + f_1(t, y_1[n]) * dt**2/2
        # solution to the 3rd order
        y_2[n+1] = y_2[n] + f_0(t, y_2[n]) * dt + f_1(t, y_2[n]) * dt**2/2 \
                          + f_2(t, y_2[n]) * dt**3/6

    # correct the final value if wanted
    if exact_end is True:
        for y in [y_0, y_1, y_2]:
            grad = (y[-1] - y[-2]) / dt
            dx = t_max - times[-2]
            y_final = y[-2] + dx * grad  # y value at t_max

            y[-1] = y_final
            times[-1] = t_max

    return times, y_0, y_1, y_2


def partial_taylor(f, f_y, f_yy, f_t, f_tt, f_yt,
                   y0, t0, t_max, dt, exact_end=False):
    """
    here you pass in the partial derivaties and through use of the
    multivairable taylor expansion formula the total derivatives are
    calculated
    """

    # define total derviative equations
    def f_1(t, y):
        """second derivative from the derivative function f and its
        partial derivatives f_y and f_t"""
        return f_t(t, y) + f(t, y) * f_y(t, y)

    def f_2(t, y):
        """first derivative from the derivative function f and its
        partial derivatives f_y, f_yy, f_t, f_tt, f_yt"""
        tmp1 = f(t, y)
        tmp2 = f_y(t, y)
        return f_tt(t, y) + 2 * tmp1 * f_yt(t, y) + f_t(t, y) * tmp2 + \
            tmp1 * tmp2**2 + tmp1**2 * f_yy(t, y)

    # pass total derivatives into the solving function
    return simple_taylor(f, f_1, f_2, y0, t0, t_max, dt, exact_end)


# Task 1C


def RK4_alpha(f, y0, t0, t_max, dt, alpha=0.5):
    """
    f is ordered f(t, y)

    aka trapezodial rule

    y0 must be list like

    exact_end : will interpolate the last point linearly to evaluate the
    funciton at exactly t_max, warning the final interval will be less than
    dt in size.
    """

    # validate inputs
    assert callable(f), 'f must be a function'
    assert alpha > 0, 'alpha must be positive'
    y_init = np.array(y0)
    assert y_init.shape != (), 'y0 must be iterable (list, array or tuple)'
    t0, t_max, dt, alpha = float(t0), float(t_max), float(dt), float(alpha)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0] = y_init

    # solve for every t
    for n, t in enumerate(times[:-1]):
        k1 = dt * f(t, y[n])
        k2 = dt * f(t + alpha*dt, y[n] + alpha*k1)
        k3 = dt * f(t + alpha*dt, y[n] + alpha*k2)
        k4 = dt * f(t + dt, y[n] + k3)
        y[n+1] = y[n] + (k1/6 + k2/3 + k3/3 + k4/6)

    return times, y


def AB4(f, y0, t0, t_max, dt):

    # validate inputs
    assert callable(f), 'f must be a function'
    y_init = np.array(y0)
    assert y_init.shape != (), 'y0 must be iterable (list, array or tuple)'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0:4] = RK4_alpha(f, y0, t0, 3*dt, dt)[1]  # initialise first 3 steps
    f_store = [f(times[0], y[0]), f(times[1], y[1]), f(times[2], y[2])]

    # solve for every t
    for n, t in enumerate(times[:-1]):
        # skip the iterations done by rk4
        if n < 3:
            continue
        f_n = f(t, y[n])
        # AB4 step
        tmp = (-9 * f_store[0] + 37 * f_store[1] - 59 * f_store[2] +
               55 * f_n) / 24
        y[n+1] = y[n] + dt * tmp
        # update f_store
        f_store = f_store[1:] + [f_n]

    return times, y


def AM3_PC(f, y0, t0, t_max, dt):

    # validate inputs
    assert callable(f), 'f must be a function'
    y_init = np.array(y0)
    assert y_init.shape != (), 'y0 must be iterable (list, array or tuple)'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0:4] = RK4_alpha(f, y0, t0, 3*dt, dt)[1]  # initialise first 3 steps
    f_store = [f(times[0], y[0]), f(times[1], y[1]), f(times[2], y[2])]

    # solve for every t
    for n, t in enumerate(times[:-1]):
        # skip the iterations done by rk4
        if n < 3:
            continue
        f_n = f(t, y[n])
        # AB4 step
        y_tmp = y[n] + dt * (-9 * f_store[0] + 37 * f_store[1] - 59 *
                             f_store[2] + 55 * f_n) / 24
        # AM3 step
        tmp2 = (f_store[1] - 5 * f_store[2] + 19 * f_n
                + 9 * f(t + dt, y_tmp)) / 24
        y[n+1] = y[n] + dt * tmp2
        # update f_store
        f_store = f_store[1:] + [f_n]

    return times, y


# Task 1D


def Euler_Backward(f, y0, t0, t_max, dt, atol=0.0001):
    "fixed point"

    # validate inputs
    assert callable(f), 'f must be a function'
    y_init = np.array(y0)
    assert y0.shape != (), 'y0 must be iterable'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0] = y_init

    # solve for every t
    for n, t in enumerate(times[:-1]):
        diff = np.array([1, 1])
        old = y[n]
        i = 0
        while all(diff > atol):
            new = y[n] + dt * f(t + dt, old)
            diff = abs(new-old)
            old = new
            i += 1
        y[n+1] = new
    return times, y


def Trapezodial(f, y0, t0, t_max, dt, atol=0.0001):
    "fixed point"

    # validate inputs
    assert callable(f), 'f must be a function'
    y_init = np.array(y0)
    assert y0.shape != (), 'y0 must be iterable'
    t0, t_max, dt = float(t0), float(t_max), float(dt)

    # set up solution
    times = np.arange(t0, t_max + dt, dt)  # 1 step over t_max
    y = np.zeros([times.size, y_init.size])  # to hold solutions in
    y[0] = y_init

    # solve for every t
    for n, t in enumerate(times[:-1]):
        # backwards step
        diff = np.array([1, 1])
        old = y[n]
        i = 0
        deriv = f(t, y[n])
        while all(diff > atol):
            new = y[n] + dt * (f(t + dt, old) + deriv) / 2
            diff = abs(new-old)
            old = new
            i += 1
        # forwards step
        y[n+1] = new
    return times, y
