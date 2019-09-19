# Numerical_Methods

A collection of numerical methods projects I have gathered through course assignment in my degrees. Each directory has more information inside on the project aims and results. Most code is implemented in python or C++.

## Repository Structure

- Airfoil_optimisation - Previous coursework example on optimising airfoil shape to have an optimal lift area ratio. Implemented in python through a Jupyter notebook and accompanying files.

- Climate_change_model - Using linear regression to predict the average monthly temperature in Tokyo. Implemented in C++ with graphical outputs from GNU plot.

- FD_FEM_methods - A small collection of numerical integrators (e.g. RK2) and practice of solving the advection-diffusion equation in either steady or time-varying states. All methods and their solutions are explained by an accompanying notebook.

- Interpolation - Holds a large collection of numerical exercises including interpolation, LU decomposition (and other matrix methods), fast fourier transforms, investigation of pseudo-random numbers and machine precision. These are all given as flexible python functions.

- Matrix_solver - A positive definite matrix solver implemented in C++ for fast execution. Memory and data structures were considered for optimal performance.

- Medical_imaging - Small utility for analysing medical (DICOM) images in C++. Has a command-line interface that allows for edge detection, colour inversion and several other filter options.

- Meteorite_simulator - Collection of python numerical methods used to simulate the altitude at which a meteorite will air burst (disintegrate into a shock wave). Code features a kwarg front end for easy and adaptable simulation, including changing the gravity and atmosphere to say simulate on mars.

- Method_of_Manufactured_Solutions - Contains an example notebook whereby the method of manufactured solutions is used to verify a numerical solver for the advection-diffusion problem. There is also a notebook exploring the use of sympy which is very useful for creating a manufactured solution.

- Numerical_Integration - Contains monte carlo and newton coats methods for integration of functions. The vegas (adaptive monte carlo method) has also been implemented and extensively investigated to show good performance on high dimensional functions.
