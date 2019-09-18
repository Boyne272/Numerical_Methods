# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 15:36:21 2018

@author: yeliu
"""
import numpy as np
import unittest as ut
from asteroid_solver import asteroid


class Test_asteroid_solver(ut.TestCase):

    def test_simulation_setup(self):
        "test simulation correctly sets up the system"
        conditions = [19/2, 19.2e3, 3300, np.inf, 18.3*np.pi/180]
        system = asteroid(*conditions)

        assert hasattr(system, 'r'), 'system incorrectly setup'
        assert hasattr(system, 'v'), 'system incorrectly setup'
        assert hasattr(system, 'dens'), 'system incorrectly setup'
        assert hasattr(system, 'stren'), 'system incorrectly setup'
        assert hasattr(system, 'angle'), 'system incorrectly setup'
        assert hasattr(system, 'height'), 'system incorrectly setup'
        assert hasattr(system, 'mass'), 'system incorrectly setup'
        assert hasattr(system, 'rhoa'), 'system incorrectly setup'
        assert hasattr(system, 'Cd'), 'system incorrectly setup'
        assert hasattr(system, 'Ch'), 'system incorrectly setup'
        assert hasattr(system, 'Cl'), 'system incorrectly setup'
        assert hasattr(system, 'alpha'), 'system incorrectly setup'
        assert hasattr(system, 'g'), 'system incorrectly setup'
        assert hasattr(system, 'Rp'), 'system incorrectly setup'

    def test_full_chelyabinsk(self, tolerance=0.1):
        """
        run the full chelyabinsk simulation
        Tolerance is not as precise here as simulation is not an exact
        match to reality
        """
        conditions = [8.5, 19.2e3, 3300, 4e6, 18.3*np.pi/180]
        system = asteroid(*conditions)
        system.solve_ode()
        peak = [system.burst['height (km)'], system.burst['ke Lost (kt/km)']]
        peak = np.array(peak)

        obs_peak = np.array([29.5578, 81.505])
        diff = abs(peak - obs_peak)
        assert np.all(diff < obs_peak * tolerance), \
            "chelyabinsk simulation does not match imperical data"

    def simplified_chelyabinsk(self, solver='ode', rtol=None):

        # setup the problem
        # conditions = [radius, velocity, density, strength, angle (rads)]
        conditions = [19/2, 19.2e3, 3300, np.inf, 18.3*np.pi/180]
        system = asteroid(*conditions)
        system.g = 0.
        system.Cl = 0.
        system.Ch = 0.
        system.Rp = np.inf

        # pick the solver
        if solver == 'ode':
            system.solve_ode(rtol=rtol)
        elif solver == 'ivp':
            system.solve_ivp(rtol=rtol)
        else:
            raise ValueError('solver name not recognised')

        return system.result[:, 0], system.result[:, 3]

    def analytic_chelyabinsk(self, z, y0):
        """
        get the analytic solution to the problem,
        y0 are the initial conditions
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

    def test_analytic_simplifed(self, tolerance=0.001):
        """
        Numerical comparison of odeint to simplifed case
        Tolerence is the relative difference allowed between the numeric
        and analytic solution
        """

        # solve for the different cases
        nv1, nz1 = self.simplified_chelyabinsk(solver='ode', rtol=1e-5)
        av = self.analytic_chelyabinsk(nz1, [19.2e3, 95e3])

        diff = abs(nv1 - av)
        assert np.all(diff < tolerance * av), \
            'Odeint does not agree with the analytic case'


if __name__ == "__main__":
    ut.main()
