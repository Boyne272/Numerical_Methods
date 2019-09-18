# -*- coding: utf-8 -*-
"""
A module built for simulating asteroids passing though the atmosphere to
assess the hazards associated with meteorites.

Refer to Doc_string or acse4-project1-report.ipynb for user manual

Created on Mon Nov  5 10:59:46 2018
Authors: Richard Boyne, Chirayu Khimji, Ye Liu, Chen Zongpeng
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci
from scipy.interpolate import interp1d


class asteroid():
    """
    Define a class called asteroid which has basic parameters for asteroid as
    input in init function and some useful functions to compute results
    """

    def __init__(self, radius, entry_velocity, density, metior_bulk_strength,
                 angle, atmos_fun=None, init_height=95e3, init_mass=None):
        """
        Sets up initial parameters and defualt constants for the system

        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters

        entry_velocity : float
            The entery speed of the asteroid in meters/second

        density : float
            The density of the asteroid in kg/m^3

        metior_bulk_strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2

        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            in radians

        Optional
        --------
        atmos_fun : callable(z)
            Computes atmospheric density, rho, at height, z
            Defualt is rho = 1.2 exp(-z/8000)

        init_height : float
            The height of the asteroid at the given initial conditions
            Default is 95km, the approximate edge of Earths atmosphere

        init_mass : float
            The intial mass of the asteroid
            Allows for more complex geometries to be modelled, however
            other assumptions are likely false in such a case
            Defualt asteroid is assumed a sphere with a uniform density that
            is given above

        Returns
        -------
        None
        """

        # Set The Passed Init Conditions
        self.r = radius
        self.v = entry_velocity
        self.dens = density
        self.stren = metior_bulk_strength
        self.angle = angle
        self.height = init_height

        # initialse later set vairables
        self.result, self.burst, self.impact = None, None, None

        # set the mass
        if init_mass is None:
            self.mass = np.pi * self.r**3 * self.dens * 4 / 3
        else:
            self.mass = init_mass

        # set atmoshper
        if callable(atmos_fun):
            self.rhoa = atmos_fun
        else:
            self.rhoa = lambda z: 1.2 * np.exp(-z/8000)

        # set constants
        self.Cd = 1.
        self.Ch = 0.1
        self.Q = 1e7
        self.Cl = 1e-3
        self.alpha = 0.3
        self.Rp = 6371e3
        self.g = 9.81

        # mis.
        self.kt_cov = 1 / 4.184e12

    def outcome(self):
        """
        Displays the Airburst and Impact data of an already solved system

        Requiers either self.solve_ode() or self.solve_ivp() to be run first

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # check a solve function has been run
        assert np.any(self.result) is not None, \
            'Solve method has not been called yet'

        # Airburst outcome
        print('\n\nAirburst')
        for key in self.burst:
            print('\n   %s :' % key, self.burst[key])

        # Createring outcome
        print('\n\nImpact details')
        for key in self.impact:
            print('\n   %s :' % key, self.impact[key])

    def solve_ode(self, n=10000, t_final=500, rtol=None, atol=None):
        """
        Solve the system with the scipy.integrate.odeint() function

        Optional Parameters
        ----------
        n : int
            The number of time values to store the vairable values
            This is not the time step resolution as a vairable dt
            is chosen by odeint

        t_final : float
            The time to solve until
            If too short a time is chosen then the asteroid may never
            reach the ground

        rtol, atol : float
            The relative and absolute tolerence to be passed to odeint
            (see scipy.integrate docstring for more details)

        Returns
        -------
        Result : 2d array of size (n, number_of_vairables)
            Solution to the system in order:
            Result = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
            (R = v, m, theta, z, x, r)
            Can also be accessed as self.result
        """

        # setup and solve
        self.times = np.linspace(0., t_final, n)
        y0 = np.array([self.v, self.mass, self.angle,
                       self.height, 0., self.r])
        self.result = sci.odeint(self.equations, y0, self.times,
                                 tfirst=True)

        # solve with rtol and atol if given
        if rtol is None and atol is None:
            self.result = sci.odeint(self.equations, y0, self.times,
                                     tfirst=True)
        elif rtol is not None and atol is None:
            self.result = sci.odeint(self.equations, y0, self.times,
                                     tfirst=True, rtol=rtol)
        elif rtol is None and atol is not None:
            self.result = sci.odeint(self.equations, y0, self.times,
                                     tfirst=True, atol=atol)
        else:
            self.result = sci.odeint(self.equations, y0, self.times,
                                     tfirst=True, rtol=rtol, atol=atol)

        # find the outcome
        self.analyse()
        return self.result

    def solve_ivp(self, n=10000, t_final=500, rtol=None, atol=None,
                  solver='RK45'):
        """
        Solve the system with the scipy.integrate.solve_ivp() function

        Optional Parameters
        ----------
        n : int
            The number of time values to store the vairable values
            This is not the time step resolution as a vairable dt
            is chosen by odeint

        t_final : float
            The time to solve until
            If too short a time is chosen then the asteroid may never
            reach the ground

        rtol, atol : float
            The relative and absolute tolerence to be passed to odeint
            (see scipy.integrate docstring for more details)

        solver : string
            The solver for scipy.integrate.solve_ivp to use (see scipy
            docstring for more details)
            Defualt is 'RK45'

        Returns
        -------
        Result : 2d array of size (n, number_of_vairables)
            Solution to the system in order:
            Result = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
            (R = v, m, theta, z, x, r)
            Can also be accessed as self.result
        """

        # setup and solve
        self.times = np.linspace(0., t_final, n)
        y0 = np.array([self.v, self.mass, self.angle,
                       self.height, 0., self.r])

        # solve with rtol and atol if given
        if rtol is None and atol is None:
            self.result = sci.solve_ivp(self.equations, (0., t_final), y0,
                                        method=solver, t_eval=self.times)
        elif rtol is not None and atol is None:
            self.result = sci.solve_ivp(self.equations, (0., t_final), y0,
                                        method=solver, t_eval=self.times,
                                        rtol=rtol)
        elif rtol is None and atol is not None:
            self.result = sci.solve_ivp(self.equations, (0., t_final), y0,
                                        method=solver, t_eval=self.times,
                                        atol=atol)
        else:
            self.result = sci.solve_ivp(self.equations, (0., t_final), y0,
                                        method=solver, t_eval=self.times,
                                        rtol=rtol, atol=atol)
        self.result = self.result['y'].T

        # find the outcomes and print them
        self.analyse()
        return self.result

    def equations(self, t, vals):
        """
        The system of differential equations to be solved

        Parameters
        ----------
        t : float
            The current time value
            (not relevant for this system of ODEs as the system is time
            invariant i.e. if t0=0 or t0=10 the outcome is the same)

        Vals : ndarray
            Array of all vairables at the current time step in order:
            vals = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
            (vals = v, m, theta, z, x, r)

        Returns
        -------
            Detivatives : array
                The time detivatives at the current vairable values
        """

        # define local names and repeatedly used values
        v, m, theta, z, x, r = vals
        A = np.pi * r**2
        a_dens = self.rhoa(z)

        # find the detrivaties
        dv = -(self.Cd * a_dens * A * v**2 / (2 * m)) \
            + self.g * np.sin(theta)
        dm = -self.Ch * a_dens * A * v**3 / (2 * self.Q)
        d_theta = (self.g * np.cos(theta) / v) \
            - self.Cl * a_dens * A * v / (2 * m) \
            - v * np.cos(theta) / (self.Rp + z)
        dz = -v * np.sin(theta)
        dx = v * np.cos(theta) / (1 + (z / self.Rp))

#        # radius shrinkage for spherical metior
#        dr = dm / (4 * np.pi * r**2 * self.dens)

        # resolve fragmentation
        if a_dens * v**2 > self.stren:
            dr = v * np.sqrt(self.alpha * 7 * a_dens / (2 * self.dens))
        else:
            dr = 0

        return np.array([dv, dm, d_theta, dz, dx, dr])

    def plot(self, opt, ax=None):
        """
        Plot the solution to the system

        Parameters
        ----------

        opt : string
            Option for which vairables to plot:

            ===========   ================================================
            opt           description
            ===========   ================================================
            'energy'      Energy vs Height
            'eph'         Energy lost per unit length vs Height
            'speed'       Velocity vs Height
            'speed-time'  Time vs Velocity
            'mass'        Mass vs Height
            'radius'      Radius vs Height
            'angle'       Trajecttory angle vs Height
            'height'      Time vs Height
            'traject'     Horizontal distance vs Height
            'all'         Plots 'energy', 'speed', 'mass' and 'radius'
                          Requiers ax to given as a (2, 2) array of axes

        ax : matplotlib axis object or array of axes
            The axis to plot on

        Returns
        -------
        None
        """

        # check a solve function has been run
        assert np.any(self.result) is not None, \
            'Solve method has not been called yet'

        # if no axes passed make one
        if np.any(ax) is None:
            fig, ax = plt.subplots()

        if opt == 'energy':
            ax.plot(self.energy * self.kt_cov, self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'Energy [kt]',
                   ylabel='Height [km]')

        elif opt == 'eph':
            ax.plot(self.eph * self.kt_cov * 1e3, self.avg_height * 1e-3)
            ax.set(xlabel=r'Energy per Unit Height [kt km$^{-1}$]',
                   ylabel='Height [km]')

        elif opt == 'speed':
            ax.plot(self.result[:, 0] * 1e-3, self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'velocity [km s$^{-1}$]',
                   ylabel='Height [km]')

        elif opt == 'speed-time':
            ax.plot(self.times, self.result[:, 0] * 1e-3)
            ax.set(xlabel='time [s]',
                   ylabel=r'velocity [km s$^{-1}$]')

        elif opt == 'mass':
            ax.plot(self.result[:, 1], self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'Mass [kg]',
                   ylabel='Height [km]')

        elif opt == 'radius':
            ax.plot(self.result[:, 5], self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'Radius [m]',
                   ylabel='Height [km]')

        elif opt == 'angle':
            ax.plot(self.result[:, 2], self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'angle [rad]',
                   ylabel='Height [km]')

        elif opt == 'height':
            ax.plot(self.times, self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'Time [s]',
                   ylabel='Height [km]')

        elif opt == 'traject':
            ax.plot(self.result[:, 4] * 1e-3, self.result[:, 3] * 1e-3)
            ax.set(xlabel=r'Horizontal Distance [km]',
                   ylabel='Height [km]')
            ax.grid()

        if opt == 'all':
            assert isinstance(ax, type(np.array([]))), \
                   'ax must be a (2, 2) array of axes'
            self.plot('energy', ax[0][0])
            self.plot('speed', ax[0][1])
            self.plot('mass', ax[1][0])
            self.plot('radius', ax[1][1])
            for axs in ax.flatten():
                axs.set(ylim=[0, self.result[0, 3]/1e3])
        else:
            # set axis limits for all single axes
            ax.set_ylim(bottom=0)
        return

    def analyse(self):
        """
        Inspect a prefound solution to calculate the impact and blast stats
        These results are stored in self.blast and self.impact, and can be
        viewed with self.outcome()

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # check a solve function has been run
        assert np.any(self.result) is not None, \
            'Solve method has not been called yet'

        # truncate results for negative altitude
        try:
            hit_ground = np.where(self.result[:, 3] < 0)[0][0] - 1
            self.result = self.result[:hit_ground]
            self.times = self.times[:hit_ground]
        except TypeError:
            print('WARNING asteroid never reached the ground, simulation \
                  incomplete. Try solving for a longer time')
        except IndexError:
            print('WARNING asteroid never reached the ground, simulation \
                  incomplete. Try solving for a longer time')

        # find the energy at each time step
        self.energy = 0.5 * self.result[:, 1] * self.result[:, 0]**2

        # finite difference to find dE/dt before hitting the ground
        energy_diff = self.energy[:-1] - self.energy[1:]
        hight_diff = self.result[:-1, 3] - self.result[1:, 3]
        self.eph = energy_diff / hight_diff
        self.avg_height = (self.result[:-1, 3] + self.result[1:, 3]) / 2

        # find the airburst height
        i_burst = self.eph.argmax()
        ke_diff = self.energy[0] - self.energy[i_burst]
        self.burst = {'occured': self.avg_height[i_burst] > 0,
                      'height (km)': self.avg_height[i_burst] / 1e3,
                      'ke Lost (kt/km)': self.eph[i_burst] * self.kt_cov * 1e3,
                      'ke Lost by burst (kt/km)': ke_diff * self.kt_cov * 1e3,
                      '% ke Lost by burst': 100*ke_diff/self.energy[0]}

        # find impact stats if ground was hit
        if 'hit_ground' in locals():
            i_impact = abs(self.result[:, 3]).argmin()
            self.impact = {'ground reached': True,
                           'time (s)': self.times[i_impact],
                           'mass (kg)': self.result[i_impact, 1],
                           'speed (km/s)': self.result[i_impact, 0] / 1e3}
        else:
            self.impact = {'ground reached': False}


def solve_ensumble(ppfs, n=1000, seed=None):
    """
    Run asteroid simulation for a distribution of initial conditions and
    find the burst distribution

    Parameters
    ----------
    ppfs : array like of callable(x)
        List of ppf functions for each initial condition (inverted cdf
        functions). These functions must take a single argument 0 <= x <= 1
        (random variable) and return a value for the initial confifition
        on the desired distribution. Must be length 5 and have
        order that follows the input for asteroid class
        (i.e. radius, entry_velocity, density, metior_bulk_strength, angle)

    Optional
    --------
    n : int
        Number of times to sample the given distributions

    seed : float
        Seed for the random number generator

    Returns
    -------

    distribution : tubple of arrays
        First array is the average for burst altitude and magnitude
        Second array is the uncertainty in each repectivly
    """

    # set the seed
    np.random.seed(seed)

    alt = []
    mag = []

    # solve for distrubuted initial conditions
    for i in range(n):
        params = [ppf(np.random.rand()) for ppf in ppfs]
        tmp_sys = asteroid(*params)
        tmp_sys.solve_ode()
        alt.append(tmp_sys.burst['height (km)'])
        mag.append(tmp_sys.burst['ke Lost (kt/km)'])

    # find the averages and vairations
    alt_un = np.std(alt)
    alt_avg = np.mean(alt)
    mag_un = np.std(mag)
    mag_avg = np.mean(mag)

    print('peak height: %.7f +- %.7f' % (alt_avg, alt_un))
    print('peak magnitude: %.7f +- %.7f' % (mag_avg, mag_un))
    return np.array([alt_avg, mag_avg]), np.array([alt_un, mag_un])


def load_atmos(file, plot=True):
    """
    Load atmospheric data from a table, linearly interpolat it and create an
    atmospheric density function

    Parameters
    ----------
    file : string
        Relative file path to csv file for atmospheric density,
        Table needs to be of the shape
            | height  |  rho0  |  H |
        where rho0 and H are parameters in the following function:
            density = rho0 * exp(-z/H)

    Optional
    --------
    plot : boolian
        Option to plot the new density function, compared to earths simple
        atmosphere model. Good to check data loaded correctly

    Returns
    -------
    f_dens : callable(z)
        returns the density at that altitude z

    where rho0 and H are parameters in the following function
        density = rho0 * exp(-z/H)

    plot is to see what the loaded atmospheric density function looks
    like incase there is an error in the way it is loaded
    """

    # load the file and interpolate
    para_list = np.loadtxt(file, delimiter=' ', unpack=True)

    inter_rho0 = interp1d(para_list[0], para_list[1])
    inter_H = interp1d(para_list[0], para_list[2])

    def f(z):
        # deal with negative heights
        if np.any(z < 0):
            return inter_rho0(0)
        rho0 = inter_rho0(z)
        H = inter_H(z)
        return rho0 * np.exp(-z/H)

    # plot new density function and simple atmoshpere
    if plot is True:
        # set up figures
        fig = plt.figure(figsize=(15, 10))
        fig.suptitle('Interpolated Atmos density Table')
        ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)
        ax2 = plt.subplot2grid((2, 2), (0, 1))
        ax3 = plt.subplot2grid((2, 2), (1, 1))

        z_vals = np.linspace(para_list[0, 0], para_list[0, -1], 100000)

        # plot found function
        ax1.plot(f(z_vals), z_vals, label='Tabulated')
        ax1.set(xlabel=r'Atmospheric Density (kg/m$^3$)', ylabel='height (m)')
        simple = 1.2*np.exp(-z_vals/8000)
        ax1.plot(simple, z_vals, '--', label='Simple model')
        ax1.legend()

        # plot read H and rho values
        ax2.plot(inter_H(z_vals), z_vals)
        ax2.set(xlabel='H (m)', ylabel='height (m)')
        ax3.plot(inter_rho0(z_vals), z_vals)
        ax3.set(xlabel=r'$\rho_0$ (kg/m$^3$)', ylabel='height (m)')

    return f


if __name__ == '__main__':
    # set matplotlib preferences - viewing
    plt.rc('axes', titlesize=20, labelsize=20)
    plt.rc('axes.formatter', limits=[-4, 4])
    plt.rc('ytick', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('lines', linewidth=1.5, markersize=7)
    plt.rc('figure', figsize=(9, 9))
    plt.rc('legend', fontsize=15)
    plt.close('all')

    # example solution
    conditions = [8.5, 19.2e3, 3300, 4.0e6, 18.3*np.pi/180]
#   conditions = [radius, velocity, density, strength, angle]
    system = asteroid(*conditions)

    # solve and plot
    system.solve_ode(rtol=1e-12)
    fig, ax1 = plt.subplots(2, 2, figsize=[12, 12])
    system.plot('all', ax1)
    fig, ax2 = plt.subplots()
    system.plot('eph', ax2)

    # plot Chelyabinsk
    cher_raw = np.loadtxt('data/ChelyabinskEnergyAltitude.csv',
                          dtype=type(''), delimiter=',')
    cher = cher_raw[1:].astype(float)
    ax2.plot(cher[:, 1], cher[:, 0], label='Measured Data')
