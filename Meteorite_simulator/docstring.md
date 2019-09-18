
Docstring for module asteriod_solver

NAME
    asteroid_solver

DESCRIPTION
    A module built for simulating asteroids passing though the atmosphere to
    assess the hazards associated with meteorites.

    Refer to Doc_string or acse4-project1-report.ipynb for user manual

    Created on Mon Nov  5 10:59:46 2018
    Authors: Richard Boyne, Chirayu Khimji, Ye Liu, Chen Zongpeng

CLASSES
    builtins.object
        asteroid

    class asteroid(builtins.object)
     |  Define a class called asteroid which has basic parameters for
     |  asteroid as input in init function and some useful functions
     |  to compute results
     |
     |  Methods defined here:
     |
     |  __init__(self, radius, entry_velocity, density, metior_bulk_strength, angle, atmos_fun=None, init_height=95000.0, init_mass=None)
     |      Sets up initial parameters and defualt constants for the system
     |
     |      Parameters
     |      ----------
     |      radius : float
     |          The radius of the asteroid in meters
     |
     |      entry_velocity : float
     |          The entery speed of the asteroid in meters/second
     |
     |      density : float
     |          The density of the asteroid in kg/m^3
     |
     |      metior_bulk_strength : float
     |          The strength of the asteroid (i.e. the maximum pressure it can
     |          take before fragmenting) in N/m^2
     |
     |      angle : float
     |          The initial trajectory angle of the asteroid to the horizontal
     |          in radians
     |
     |      Optional
     |      --------
     |      atmos_fun : callable(z)
     |          Computes atmospheric density, rho, at height, z
     |          Defualt is rho = 1.2 exp(-z/8000)
     |
     |      init_height : float
     |          The height of the asteroid at the given initial conditions
     |          Default is 95km, the approximate edge of Earths atmosphere
     |
     |      init_mass : float
     |          The intial mass of the asteroid
     |          Allows for more complex geometries to be modelled, however
     |          other assumptions are likely false in such a case
     |          Defualt asteroid is assumed a sphere with a uniform density that
     |          is given above
     |
     |      Returns
     |      -------
     |      None
     |
     |  analyse(self)
     |      Inspect a prefound solution to calculate the impact and blast stats
     |      These results are stored in self.blast and self.impact, and can be
     |      viewed with self.outcome()
     |
     |      Parameters
     |      ----------
     |      None
     |
     |      Returns
     |      -------
     |      None
     |
     |  equations(self, t, vals)
     |      The system of differential equations to be solved
     |
     |      Parameters
     |      ----------
     |      t : float
     |          The current time value
     |          (not relevant for this system of ODEs as the system is time
     |          invariant i.e. if t0=0 or t0=10 the outcome is the same)
     |
     |      Vals : ndarray
     |          Array of all vairables at the current time step in order:
     |          vals = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
     |          (vals = v, m, theta, z, x, r)
     |
     |      Returns
     |      -------
     |          Detivatives : array
     |              The time detivatives at the current vairable values
     |
     |  outcome(self)
     |      Displays the Airburst and Impact data of an already solved system
     |
     |      Requiers either self.solve_ode() or self.solve_ivp() to be run first
     |
     |      Parameters
     |      ----------
     |      None
     |
     |      Returns
     |      -------
     |      None
     |
     |  plot(self, opt, ax=None)
     |      Plot the solution to the system
     |
     |      Parameters
     |      ----------
     |
     |      opt : string
     |          Option for which vairables to plot:
     |
     |          ===========   ================================================
     |          opt           description
     |          ===========   ================================================
     |          'energy'      Energy vs Height
     |          'eph'         Energy lost per unit length vs Height
     |          'speed'       Velocity vs Height
     |          'speed-time'  Time vs Velocity
     |          'mass'        Mass vs Height
     |          'radius'      Radius vs Height
     |          'angle'       Trajecttory angle vs Height
     |          'height'      Time vs Height
     |          'traject'     Horizontal distance vs Height
     |          'all'         Plots 'energy', 'speed', 'mass' and 'radius'
     |                        Requiers ax to given as a (2, 2) array of axes
     |
     |      ax : matplotlib axis object or array of axes
     |          The axis to plot on
     |
     |      Returns
     |      -------
     |      None
     |
     |  solve_ivp(self, n=10000, t_final=500, rtol=None, atol=None, solver='RK45')
     |      Solve the system with the scipy.integrate.solve_ivp() function
     |
     |      Optional Parameters
     |      ----------
     |      n : int
     |          The number of time values to store the vairable values
     |          This is not the time step resolution as a vairable dt
     |          is chosen by odeint
     |
     |      t_final : float
     |          The time to solve until
     |          If too short a time is chosen then the asteroid may never
     |          reach the ground
     |
     |      rtol, atol : float
     |          The relative and absolute tolerence to be passed to odeint
     |          (see scipy.integrate docstring for more details)
     |
     |      solver : string
     |          The solver for scipy.integrate.solve_ivp to use (see scipy
     |          docstring for more details)
     |          Defualt is 'RK45'
     |
     |      Returns
     |      -------
     |      Result : 2d array of size (n, number_of_vairables)
     |          Solution to the system in order:
     |          Result = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
     |          (R = v, m, theta, z, x, r)
     |          Can also be accessed as self.result
     |
     |  solve_ode(self, n=10000, t_final=500, rtol=None, atol=None)
     |      Solve the system with the scipy.integrate.odeint() function
     |
     |      Optional Parameters
     |      ----------
     |      n : int
     |          The number of time values to store the vairable values
     |          This is not the time step resolution as a vairable dt
     |          is chosen by odeint
     |
     |      t_final : float
     |          The time to solve until
     |          If too short a time is chosen then the asteroid may never
     |          reach the ground
     |
     |
     |      rtol, atol : float
     |          The relative and absolute tolerence to be passed to odeint
     |          (see scipy.integrate docstring for more details)
     |
     |      Returns
     |      -------
     |      Result : 2d array of size (n, number_of_vairables)
     |          Solution to the system in order:
     |          Result = Velocity, Mass, Angle, Height, Horizontal Distance, Radius
     |          (R = v, m, theta, z, x, r)
     |          Can also be accessed as self.result
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

FUNCTIONS
    load_atmos(file, plot=True)
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

    solve_ensumble(ppfs, n=1000, seed=None)
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