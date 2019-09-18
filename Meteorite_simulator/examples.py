#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
examples of test_asteroid_solver at work

Created on Thu Nov  8 09:23:11 2018
Authors: Richard Boyne, Chirayu Khimji, Ye Liu, Chen Zongpeng
"""

from asteroid_solver import asteroid, load_atmos, solve_ensumble
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.stats as scip


def chelyabinsk_demo():
    """
    Simulate the chelyabinsk passing through the atmosphere
    """

    conditions = [19/2, 19.3e3, 3300, 1e6, 18.3*np.pi/180]
    system = asteroid(*conditions)
    system.solve_ode()

    fig, ax = plt.subplots(1, 1, figsize=[12, 12])
    system.plot('eph', ax)
    ax.lines[-1].set_label('numerical solution')
    ax.set_title('Chelyabinsk Energy Differential')

    # plot Chelyabinsk
    cher_raw = np.loadtxt('data/ChelyabinskEnergyAltitude.csv',
                          dtype=type(''), delimiter=',')
    cher = cher_raw[1:].astype(float)
    ax.plot(cher[:, 1], cher[:, 0], label='Measured Data')


def simplified_chelyabinsk(rtol, solver, ret='vz'):
    """
    Solve for v and z in the simplified chelyabinsk case with
    a choice of solver ('ode' or 'ivp') and rtol
    """

#    # setup the problem
#    conditions = [radius, velocity, density, strength, angle (rads)]
    conditions = [19/2, 19.2e3, 3300, np.inf, 18.3*np.pi/180]
    system = asteroid(*conditions)
    system.g = 0.
    system.Cl = 0.
    system.Ch = 0.
    system.Rp = np.inf

    # pick the solver and solve
    if solver == 'ode':
        system.solve_ode(rtol=rtol)
    elif solver == 'ivp':
        system.solve_ivp(rtol=rtol)
    else:
        raise ValueError('solver name not recognised')

    if ret == 'vz':
        return system.result[:, 0], system.result[:, 3]
    elif ret == 'sys':
        return system
    else:
        raise ValueError('ret value not recognised')


def analytic_chelyabinsk(z, y0):
    """
    get the analytic solution to the simplifed chelyabinsk problem,
    y0 are the initial conditions [v0, z0]
    returns v for the given z values
    """

    # set constants and ICs for chelyabinsk simplified
    Cd = 1.
    r = 19/2
    A = np.pi * r**2
    m = 3300 * np.pi * r**3 * 4 / 3
    theta = 18.3 * np.pi / 180  # in rads
    rho = 1.2
    H = 8000

    # find integration constant
    v0, z0 = y0
    tmp_cons = H * Cd * A * rho * np.exp(-z0 / H)
    C = np.log(v0) + tmp_cons / (2 * m * np.sin(theta))

    # get solution
    ln_v = C - H * Cd * A * rho * np.exp(-z / H) / (2 * m * np.sin(theta))
    return np.exp(ln_v)


def visual_comp_chelyabinsk():
    """
    visual comparison of the nueric and analytic solutions for the
    chelyabinsk case
    """

    # solve for the different cases
    nv1, nz1 = simplified_chelyabinsk(rtol=0.001, solver='ode')
    av = analytic_chelyabinsk(nz1, [19.2e3, 95e3])

    # plot the results
    fig, ax = plt.subplots()
    ax.plot(nv1/1e3, nz1/1e3, label='odeint')
    ax.plot(av[::20]/1e3, nz1[::20]/1e3, 'go', ms=5, label='analytic')
    ax.set(title='Chelyabinsk simplified numeric analtic comparison',
           xlabel='Velocity [km/s]',
           ylabel='Height [km]')
    ax.legend()
    return


def solver_comparison():
    "compare the different solvers on the non-symplified chelyabinsk case"
    conditions = [19/2, 19.2e3, 3300, 1e6, 18.3*np.pi/180]
    system = asteroid(*conditions)

    fig, ax = plt.subplots()

    system.solve_ode()
    system.plot('eph', ax)
    ax.lines[-1].set_label('odeint')

    system.solve_ivp()
    system.plot('eph', ax)
    ax.lines[-1].set_label('RK45')

    system.solve_ivp(solver='RK23')
    system.plot('eph', ax)
    ax.lines[-1].set_label('RK23')

    system.solve_ivp(solver='Radau')
    system.plot('eph', ax)
    ax.lines[-1].set_label('Radau')

    system.solve_ivp(solver='BDF')
    system.plot('eph', ax)
    ax.lines[-1].set_label('BDF')

    system.solve_ivp(solver='LSODA')
    system.plot('eph', ax)
    ax.lines[-1].set_label('LSODA')

    ax.legend()
    ax.set(xlim=[50, 80], ylim=[30, 40])


def rtol_error_simplified_chelyabinsk():
    """
    solve for a range of rtol values in the simplififed chelyabinsk
    case to gain an estimate for the error plot
    """

    # setup the wanted rtols
    res = []
    rtols = np.array([0.01 / 2**n for n in range(0, 10)])

    # solve for the different cases and find the difference
    for rtol in rtols:
        nv, nz = simplified_chelyabinsk(solver='ode', rtol=rtol)
        av = analytic_chelyabinsk(nz, [19.2e3, 95e3])
        diff = abs(nv-av).mean()
        res.append(diff)

    # find the error order
    x, y = np.log10(rtols), np.log10(np.array(res))
    params = np.polyfit(x, y, 1)
    print(params)

    # plot the result
    fig, ax = plt.subplots()
    ax.loglog(rtols, res, 'bo', label='Error')
    ax.loglog(rtols, 10**params[1] * rtols**params[0], 'r--', label='Order')
    ax.set(title='Avg Error vs Rtol for Chelyabinsk Simplified',
           xlabel='Rtol', ylabel='Avg Error [m]')
    ax.legend()


def find_chelyabinsk_params():
    """
    find the parameters asteroid size and strength
    This is done by minimising the difference from the measured and calculated
    values for Chelyabinsk using the data in
    """

    # load chelyabinsk data and find the peak
    cher_raw = np.loadtxt('data/ChelyabinskEnergyAltitude.csv', dtype=type(''),
                          delimiter=',')
    cher = cher_raw[1:].astype(float)
    peak = cher[cher[:, 1].argmax()]

    # setup the system
    conditions = [8.2, 19.2e3, 3300, 4e6, 18.3*np.pi/180]
    system = asteroid(*conditions)
    init_guess = np.array([conditions[0], conditions[3]])

    def to_minimise(params, system, peak):
        # set the new params
        system.radius, system.stren = params
        system.mass = 4*np.pi*system.radius**3 * system.dens / 3

        # solve the system
        system.solve_ode()
        system.analyse()

        # return the difference
        diff = np.array([system.burst['height (km)'] - peak[0],
                         system.burst['ke Lost (kt/km)'] - peak[1]])
        return (diff**2).sum()

    sol = opt.minimize(to_minimise, init_guess, args=(system, peak),
                       method='Powell')

    # plot data, initial guess and optimisation
    fig, ax = plt.subplots()

#    system2 = asteroid(*conditions)
#    system2.solve_ode()
#    system2.plot('eph', ax)
#    ax.lines[-1].set(label='init guess', color='g', linestyle='--')

    conditions2 = [sol['x'][0], 19.2e3, 3300, sol['x'][1], 18.3*np.pi/180]
    system = asteroid(*conditions2)
    system.solve_ode()
    system.plot('eph', ax)
    ax.lines[-1].set(label='optimised', color='r', linestyle='-')

    ax.plot(cher[:, 1], cher[:, 0], 'bx', ms=5,
            label='Measured Data')

    ax.set(title='Estimation of Chelyabinsk Radius & Strength',
           ylim=[20, 42])
    ax.legend()
    print('\nRadius (m) = %.3f \nStrength (Pa) = %.3f \n' %
          (sol['x'][0], sol['x'][1]))

    return sol


def compare_earth_tabulated():

    # setup the systems
    conditions = [8.5, 19.2e3, 3300, 4e6, 18.3*np.pi/180]
    file = 'data/AltitudeDensityTable.csv'
    system_tab = asteroid(*conditions, atmos_fun=load_atmos(file),
                          init_height=8.5e4)
    system_norm = asteroid(*conditions, init_height=8.5e4)

    # solve for each
    system_tab.solve_ode()
    system_norm.solve_ode()

    # plot the trajectories
    fig, ax_tra = plt.subplots()
    system_tab.plot('traject', ax_tra)
    ax_tra.lines[-1].set_label('Tabulated Atmosphere')
    system_norm.plot('traject', ax_tra)
    ax_tra.lines[-1].set_label('Simple Atmosphere')

    # plot the energy per height
    fig, ax_eph = plt.subplots()
    system_tab.plot('eph', ax_eph)
    ax_eph.lines[-1].set_label('Tabulated Atmosphere')
    system_norm.plot('eph', ax_eph)
    ax_eph.lines[-1].set_label('Simple Atmosphere')

    # label axes
    ax_tra.set(title='Trajectory Comparison')
    ax_tra.legend()
    ax_eph.set(title='Energy Dissapation Comparison')
    ax_eph.legend()

    return system_tab


def atmdens_mars(h):
    p = 0.699*np.exp(-0.00009*h)

    if h >= 7000:
        T = -23.4 - 0.00222*h

    elif h < 7000:
        T = -31 - 0.000998*h

    return p/(0.1921*(T+273.1))


def compare_mars_tabulated():
    # setup the systems
    conditions = [8.5, 19.2e3, 3300, 4e6, 18.3*np.pi/180]
    file = 'data/AltitudeDensityTable.csv'
    system_tab = asteroid(*conditions, atmos_fun=load_atmos(file, plot=False),
                          init_height=8.5e4)
    system_mars = asteroid(*conditions, init_height=95e3,
                           atmos_fun=atmdens_mars)
    system_mars.Rp = 3390*1000
    system_mars.g = 3.8

    # solve for each
    system_tab.solve_ode()
    system_mars.solve_ode()

    # plot the trajectories
    fig, ax_tra = plt.subplots()
    system_tab.plot('traject', ax_tra)
    ax_tra.lines[-1].set_label('Tabulated Atmosphere')
    system_mars.plot('traject', ax_tra)
    ax_tra.lines[-1].set_label('Mars Atmosphere')

    # plot the energy per height
    fig, ax_eph = plt.subplots()
    system_tab.plot('eph', ax_eph)
    ax_eph.lines[-1].set_label('Tabulated Atmosphere')
    system_mars.plot('eph', ax_eph)
    ax_eph.lines[-1].set_label('Mars Atmosphere')

    # label axes
    ax_tra.set(title='Martian Trajectory vs Earth')
    ax_tra.legend()
    ax_eph.set(title='Martian Energy Dissapation vs Earth', ylim=[0, 100])
    ax_eph.legend()

    return system_tab


def ppf(n):
    """
    return the ppf velocity of a single value uniformly distributed
    from 0-1
    """
    return scip.norm.ppf(n, 19.2e3, 1000)


def test_ensumble():
    """
    test for the solve_ensumble function. Use velocity as the changing variable
    in Gaussian distribution
    the result takes the form (mean +- standard deviation)
    """
    function_list = [lambda x: 9.5, ppf, lambda x: 3300, lambda x: 4e6,
                     lambda x: 18.3*np.pi/180]
    solve_ensumble(function_list)
