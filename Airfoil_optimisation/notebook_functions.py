# these are the functions given to us in the 4 notebooks that we are not requiered to change
# they are imported from here to keep the notebook cleaner

import numpy as np
import scipy.linalg as sl
import scipy.sparse.linalg as spl
import scipy.optimize as sop
import matplotlib.pyplot as plt
import pyamg

# --------------------------- notebook 3

class Preconditioner(spl.LinearOperator):
    """A base class for the preconditioners
    
    Every preconditioner needs to implement a __call__ method
    that will be called with the residual vector and should return
    the preconditioner applied to that residual vector."""
    def __init__(self, A):
        super().__init__(dtype=None, shape=A.shape)

    def __call__(self, r):
        raise NotImplemented("The call method should be overloaded")

    def _matvec(self, r):
        return self(r)


class JacobiPreconditioner(Preconditioner):
    """Jacobi preconditioner"""
    def __init__(self, A):
        super().__init__(A)
        self.diag = A.diagonal()  # extract main diagonal

    def __call__(self, r):
        return r/self.diag


class SSORPreconditioner(Preconditioner):
    """SSOR preconditioner,
    
    set 0<omega<2 to improve convergence (see lecture notes 7)."""
    def __init__(self, A, omega=1.0):
        super().__init__(A)
        self.A = A
        self.rr = np.zeros(A.shape[0])  # auxilary vector
        self.omega = omega

    def __call__(self, r):
        self.rr[:] = 0  # reset to zero
        # this step computes rr = rr + M^{-1} r
        # since we start with rr = 0, we get rr = M^{-1} r back
        pyamg.relaxation.relaxation.sor(self.A, self.rr, r, omega=self.omega, sweep='symmetric')
        return self.rr


class AMGPreconditioner(Preconditioner):
    """AMG preconditioner - wraps around pyamg's ruge_stuben_solver"""
    def __init__(self, A):
        super().__init__(A)
        ml = pyamg.ruge_stuben_solver(A.tocsr(), strength=('classical'))
        self.gamg = ml.aspreconditioner()
        
    def __call__(self, r):
        return self.gamg(r)
    
    
class SimpleCounter:
    """Simple counter object that records how many times it has been called"""
    def __init__(self):
        self.count = 0
    
    def __call__(self, *args):
        """Counts each call, ignore any arguments"""
        self.count += 1
        
        
# --------------------------- notebook 2
        
def regular_polygon(R, n):
    """Create a regular polygon of radius R with n sides"""
    alpha = 2*np.pi/n
    x = []
    for i in range(n+1):
        x.append([R*np.cos(alpha*i), R*np.sin(alpha*i)])
    return np.array(x)


def area(x):
    """Compute the area of the polygon x"""
    n = len(x)-1
    a = 0.
    for i in range(n):
        a += (x[i,0]-x[i+1,0])*(x[i,1]+x[i+1,1])/2.
    return a


def length(x):
    """Compute the length of the path described by the list of points x."""
    n = len(x)-1
    l = 0.
    for i in range(n):
        l += sl.norm(x[i+1]-x[i])
    return l


def grad_length(x):
    """Compute the gradient of length(x)."""
    n = len(x)-1
    dldx = np.zeros((n+1, 2))
    for i in range(n):
        dx = sl.norm(x[i+1]-x[i])
        dldx[i,   0] += -(x[i+1,0]-x[i,0])/dx
        dldx[i+1, 0] +=  (x[i+1,0]-x[i,0])/dx
        dldx[i,   1] += -(x[i+1,1]-x[i,1])/dx
        dldx[i+1, 1] +=  (x[i+1,1]-x[i,1])/dx
    return dldx


# note this is actually the improved taylor test from notebook 4

def taylor_test(f, grad_f, x, h0, max_iter=20, plot_convergence=True, print_convergence=True):
    """Taylor test to verify that the function grad_f is the derivative of the function `f`
    
    We test:
    
       f(x + h) = f(x) + grad_f(x).h + O(h^2)
    
    f, grad_f  - function and its derivative to test
    x          - point in which to test
    h0         - initial perturbation
    max_iter   - number of times that h is halved
    plot_convergence - whether to plot the convergence of the Taylor residual f(x+h)-f(x)-grad_f(x).h
    print_convergence - prints the order of convergence between subsequent iterations."""


    h = h0.copy()
    residuals = []
    hnorms = []
    fx = f(x)
    gradfx = grad_f(x)
    for i in range(max_iter):
        r = f(x + h) - fx - np.tensordot(gradfx, h, len(h.shape))
        residuals.append(sl.norm(r))
        hnorms.append(sl.norm(h))
        h /= 2.

    if plot_convergence:
        fig, ax = plt.subplots(1,2, figsize=(16,4))
        ax[0].semilogy(residuals)
        ax[0].set_xlabel('iteration')
        ax[0].set_ylabel('Taylor residual $|r|$')
        ax[0].set_xticks(range(0,max_iter,2))
        ax[1].loglog(hnorms, residuals)
        ax[1].set_xlabel('perturbation size $\|h\|$')
        ax[1].set_ylabel('Taylor residual $|r|$')

        # NOTE: slope_marker does not seem to work in semilogy plots
        #annotation.slope_marker((1e-3, 1e-4), (2, 1), invert=True, ax=ax[1], size_frac=.2)

    if print_convergence:
        residuals = np.array(residuals)
        print('Order of convergence, log(r(h_i)/r(h_{i+1}))/log(2):')
        print(np.log(residuals[:-1]/residuals[1:])/np.log(2))
        
        
class ShapeOptimisationProblem:
    """Shape optimisation based on a polygon.
    
    The first and last vertex of this polygon, described a list of vertices
    should be the same and will not be changed during the optimisation."""
    def __init__(self, x0, f_user, grad_f_user, store_function_evaluations=True, store_gradient_evaluations=True):
        """
        x0     - the initial polygon of n edges, described by a shape (n+1,2) array of vertex positions,
                 where the first and last vertex should be the same to describe a closed polygon
        f_user      - the function that computes the objective values, should take a shape (n+1,2) array
        grad_f_user - function that computes the gradient of f, should take and return a shape (n+1,2) array
        store_function_evaluations, store_gradient_evaluations - if either set to true, self.x_i stores
                                                                 a list of (n+1,2) arrays representing
                                                                 the polygons for which the function or 
                                                                 gradient was evaluated."""        
        self.f_user = f_user
        self.grad_f_user = grad_f_user
        self.x0 = x0.copy()  # we keep a copy of the original polygon
        self.x = x0.copy()  # this is the copy we are actually going to change
        self.n = len(x0)-1
        self.x_i = []  # stores subsequent iterations when asked for
        self.store_function_evaluations = store_function_evaluations
        self.store_gradient_evaluations = store_gradient_evaluations
        
    def minimize(self, *args, **kwargs):
        """Wrapper around scipy.optimize.minimize
        
        Passes the correct function, gradient and initial guess. Any further arguments are also passed on
        to scipy's minimize"""
        # make a copy of x0, so that if we call minimize() again
        # we again start from the same initial guess provided in __init__()
        self.x = self.x0.copy()
        # now turn that into a flat vector
        x0_vec = self.get_x_vec()
        return sop.minimize(self.f, x0_vec, jac=self.jac, *args, **kwargs)
        
    def set_x(self, x_vec):
        """Set self.x based on the 2(n-1) vector used by scipy minimize"""
        # copy into self.x, leaving first and last vertex unchanged:
        self.x[1:self.n, :] = x_vec.reshape((self.n-1, 2))
        
    def get_x_vec(self):
        """Obtain the 2(n-1) vector from self.x"""
        return self.x[1:self.n,:].flatten()  # skip first and last vertex, and flatten
                
    def f(self, x_vec):
        """Wrapper function f, called by scipy's minimize"""
        self.set_x(x_vec)
        if self.store_function_evaluations:
            self.x_i.append(self.x.copy())
        f = self.f_user(self.x)
        return f
    
    def jac(self, x_vec):
        """Wrapper gradient function, called by scipy's minimize"""
        self.set_x(x_vec)
        if self.store_gradient_evaluations:
            self.x_i.append(self.x.copy())
        grad = self.grad_f_user(self.x)
        return grad[1:self.n,:].flatten()  # skip gradient associated with first and last vertex
    
    
# ----------------------- These functions are written by me for question 5
#                         They are included here as functions from notebook 4 below
#                         require these, they are also written in the technical report
    
def grad_area(x):
    """Compute the gradient of area(x)."""
    n = len(x)-1
    dadx = np.zeros((n+1, 2))
    for i in range(n):
        dadx[i,   0] +=  (x[i,1] + x[i+1, 1]) / 2
        dadx[i+1, 0] += -(x[i,1] + x[i+1, 1]) / 2
        dadx[i,   1] +=  (x[i,0] - x[i+1, 0]) / 2
        dadx[i+1, 1] +=  (x[i,0] - x[i+1, 0]) / 2
    return dadx


def shape_factor(x):
    """Computes the ratio l/sqrt{A}"""
    return length(x)/np.sqrt(area(x))


def grad_shape_factor(x):
    A = area(x)
    fac1 = grad_length(x) / np.sqrt(A)
    fac2 = length(x) * grad_area(x) * (-1/2) * A**(-3/2)
    return fac1 + fac2


# ----------------------- notebook 4


class AirfoilFunctional:
    """A functional combining lift and shape_factor of an airfoil.
    
    This functional relies on the functions shape_factor(), grad_shape_factor, length(x),
    and grad_length(x) being defined elsewhere"""
    def __init__(self, verbose=True, aspect_ratio=20.0, Q_max=3.8, Q_scale=1e-2, Q_eps=1e-4):
        """
        verbose - if true print lift and shape factor and combined function value at each evaluation
        aspect_ratio - used in shape_factor: aspect ratio of ellipse that has minimal value for Q
        Q_max - maximal allowable value for shape factor Q
        Q_scale - scales the penalty parameter
        Q_eps - cap log term at this distance from Q_max
        """
        self.verbose = verbose
        self.aspect_ratio = aspect_ratio
        self.Q_max = Q_max
        self.Q_scale = Q_scale
        self.Q_eps = Q_eps
        
    def shape_factor_ar(self, x):
        """Rescaled shape_factor that is minimal for an ellipse of the specified aspect ratio."""
        x[:,1] *= self.aspect_ratio  # scale y-coordinates: ellipse -> circle
        Q = shape_factor(x)
        x[:,1] /= self.aspect_ratio   # change back y-coords
        return Q

    def grad_shape_factor_ar(self, x):
        """Gradient of shape_factor_ar"""
        x[:,1] *= self.aspect_ratio  # scale y coordinates: ellipse -> circle
        dQdx = grad_shape_factor(x)
        x[:,1] /= self.aspect_ratio  # change back y-coords
        # also rescale the y-derivatives:
        dQdx[:,1] *= self.aspect_ratio
        return dQdx
    
    def f(self, sigma, x):
        """Functional: f(sigma, x) = -lift + Q_scale * shape_factor"""
        lift = -2*sigma[-1] * length(x)
        Q = self.shape_factor_ar(x)
        penalty = - self.Q_scale * np.log(np.maximum(self.Q_max-Q, self.Q_eps))
        if self.Q_max-Q<self.Q_eps:
            penalty += self.Q_scale * (Q-self.Q_max+self.Q_eps)/self.Q_eps

        f_value = -lift + penalty
        if self.verbose:
            print("C_L, Q, penalty, f(sigma, x)", lift, Q, penalty, f_value)

        return f_value
    
    def grad_x_f(self, sigma, x):
        """Partial derivative of f with respect to x"""
        # only the shape factor term explicitly depends on x
        Q = self.shape_factor_ar(x)
        g = self.Q_scale * self.grad_shape_factor_ar(x) / np.maximum(self.Q_max-Q, self.Q_eps)
        g += 2*sigma[-1] * grad_length(x)
        return g
    
    def grad_sigma_f(self, sigma, x):
        """Partial derivative of f with respect to sigma"""
        g = np.zeros_like(sigma)
        g[-1] = 2. * length(x)
        return g
    
    
class VerticalShapeOptimisation:
    """Shape optimisation based on a polygon, in which the vertices are moved vertically only
    
    The first and last vertex of this polygon, described a list of vertices
    should be the same and will not be changed during the optimisation."""
    def __init__(self, x0, f_user, grad_f_user, store_function_evaluations=True, store_gradient_evaluations=True):
        """
        x0     - the initial polygon of n edges, described by a shape (n+1,2) array of vertex positions,
                 where the first and last vertex should be the same to describe a closed polygon
        f_user      - the function that computes the objective values, should take a shape (n+1,2) array
        grad_f_user - function that computes the gradient of f, should take and return a shape (n+1,2) array
        store_function_evaluations, store_gradient_evaluations - if either set to true, self.x_i stores
                                                                 a list of (n+1,2) arrays representing
                                                                 the polygons for which the function or 
                                                                 gradient was evaluated."""        
        self.f_user = f_user
        self.grad_f_user = grad_f_user
        self.x0 = x0.copy()  # we keep a copy of the original polygon
        self.x = x0.copy()  # this is the copy we are actually going to change
        self.n = len(x0)-1
        self.x_i = []  # stores subsequent iterations when asked for
        self.store_function_evaluations = store_function_evaluations
        self.store_gradient_evaluations = store_gradient_evaluations
        
    def minimize(self, *args, **kwargs):
        """Wrapper around scipy.optimize.minimize
        
        Passes the correct function, gradient and initial guess. Any further arguments are also passed on
        to scipy's minimize"""
        # make a copy of x0, so that if we call minimize() again
        # we again start from the same initial guess provided in __init__()
        self.x = self.x0.copy()
        # now turn that into a flat vector
        x0_vec = self.get_x_vec()
        return sop.minimize(self.f, x0_vec, jac=self.jac, *args, **kwargs)
        
    def set_x(self, x_vec):
        """Set self.x based on the (n-1)-vector used by scipy minimize"""
        # copy into self.x, leaving first and last vertex unchanged:
        self.x[1:self.n, 1] = x_vec
        
    def get_x_vec(self):
        """Obtain the (n-1)-vector from self.x"""
        return self.x[1:self.n,1].copy()
                
    def f(self, x_vec):
        """Wrapper function f, called by scipy's minimize"""
        self.set_x(x_vec)
        if self.store_function_evaluations:
            self.x_i.append(self.x.copy())
        f = self.f_user(self.x)
        return f
    
    def jac(self, x_vec):
        """Wrapper gradient function, called by scipy's minimize"""
        self.set_x(x_vec)
        if self.store_gradient_evaluations:
            self.x_i.append(self.x.copy())
        grad = self.grad_f_user(self.x)
        return grad[1:self.n,1].copy()  # skip gradient associated with first and last vertex
    
    