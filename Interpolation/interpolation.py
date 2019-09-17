import numpy as np
import matplotlib.pyplot as plt
import matrix_methods as lu

def lin_interpol(x_p, y_p):
    """
    Find the coefficents then plot the linear interpolation of (x_p, y_p) 
    points where x_p and x_p are arrays.
    """
    f = np.zeros([ x_p.shape[0] - 1 , 4 ]) # Coefficents and interval array
    
    for i in range( x_p.shape[0] - 1 ): # for every x[i], x[i+1] pair
    
        x_coeff = (y_p[i+1] - y_p[i]) / (x_p[i+1] - x_p[i])
        const = (x_p[i+1]*y_p[i] - x_p[i]*y_p[i+1] ) / (x_p[i+1] - x_p[i])
        
        # save the x coefficent, constant and the interval for this line
        f[i,:] = x_coeff, const, x_p[i], x_p[i+1]
    
    for a, b, start, end in f: # for every line fitted
        line_x = np.linspace( start, end, 3) # points to plot in x_range
        line_y = line_x * a + b # find the fitted line value at these points
        plt.plot(line_x,line_y,'k--', lw = 1, label = 'Linear' if a==f[0][0] else "") # only label one plot

def derivatives(x_p, y_p):
    """
    Finds the 2nd derivatives reuiered for the cubic interpolation equation via
    creating the matrrix equation Mx = d where x is the 2nd derivatives with M
    and d calculated from the points.
    x_p is an array of x points with corresponding y_p points.
    """
    # set up the matrix equation
    n = x_p.shape[0]
    M = np.zeros( [n,n] )
    d = np.zeros( [n,1] )
    
    # fill in the constants where they can be
    for i in np.arange(1,n-1 ): # for all but the first and last row
        M[i,i-1  ] = ( x_p[i] - x_p[i-1] ) / 6.
        M[i,i] = ( x_p[i+1] - x_p[i-1] ) / 3.
        M[i,i+1] = ( x_p[i+1] - x_p[i] ) /6.
        d[i,0  ] = ( y_p[i+1] - y_p[i] ) / ( x_p[i+1] - x_p[i] ) - ( y_p[i] - y_p[i-1] ) / ( x_p[i] - x_p[i-1] )
    
    M[0,0],M[-1,-1] = 1.,1. # compactly sets the BCs
    
    LU = lu.LU_decomp(M) # solves the matrix equations
    return lu.FB_sub(LU.Low, LU.Upp, d) # find and return 2nd derivatives

def f(X, params):
    """
    Does the cubic interpolation sum, x can be an array of points,
    params is the set of parameters for the sum, requiers:
    xi, x(1+1), yi, y(i+1), y''i, y''(i+1)
    **could be compacted but makes code too confusing**
    """
    x_i, x_i1, y_i, y_i1, y_ip, y_ip1 = params
    A = (x_i1 - X) / (x_i1 - x_i)
    B = (X - x_i) / (x_i1 - x_i)
    C = (1./6) * (A**3 - A) * (x_i1-x_i)**2 
    D = (1./6) * (B**3 - B) * (x_i1-x_i)**2
    return A*y_i + B*y_i1 + C*y_ip + D*y_ip1

def cubic_interpol(X_P, Y_P):
    """
    Plots the cubic interpolation of (X_P, Y_P) points where X_P and Y_P are
    arrays and X_P is ordered from low to high. Calls the derivatives and f functions.
    """
    y_derivs = derivatives( X_P, Y_P ).flatten() # flatten as FB_sub returns 2d array
    
    for j in np.arange( X_P.shape[0] - 1 ): # for every x[i] and x[i+1] pair
        plot_points = np.linspace( X_P[j], X_P[j+1], 20) # points to plot in the interval
        params = [ X_P[j], X_P[j+1], Y_P[j], Y_P[j+1],
                   y_derivs[j], y_derivs[j+1]]
        f_points = f(plot_points, params)
        plt.plot(plot_points, f_points, 'b-', ms = .5, label = 'Cubic'if j==0 else "") # only label one plot

def solve_i():
    """
    Solves for the range of points given both cubic and linear interpolation 
    """
    x = np.array([ -2.1, -1.45, -1.3, -0.2, 0.1, 0.15, 0.8, 1.1, 1.5, 2.8, 3.8 ])
    y = np.array([0.012155, 0.122151, 0.184520, 0.960789, 0.990050, 0.977751,
                         0.527292, 0.298197, 0.105399, 3.936690E-4, 5.355348E-7])
    # find and plot both interpolations and the oiginal points
    plt.figure(1)
    cubic_interpol(x,y)
    lin_interpol(x,y)
    plt.plot(x, y, 'rx', ms = 10, label = 'Points')
    # plot settings
    plt.title('Cubic & Linear Interpolation Given Points')
    plt.xlabel('x',fontsize = 14)
    plt.ylabel('y',fontsize = 14)
    plt.legend()

def validate_i(n = 20, r = 5, s = 0, xrand = False):
    """
    interpolate a n random y values over a range of x values seeded with s
    """
    # set the random input
    np.random.seed(s)
    x = np.arange(n)
    # randomise the x values instead of uniform spacing 
    if xrand == True:
        x = np.random.random_sample(n) * r
        x.sort()
    y = np.random.random_sample(n) * r
    # find and plot both interpolations and the oiginal points
    plt.figure(2)
    cubic_interpol(x,y)
    lin_interpol(x,y)
    plt.plot(x, y, 'rx', ms = 10, label = 'Points')
    # plot settings
    plt.title('Cubic & Linear Interpolation Random Points')
    plt.xlabel('x',fontsize = 14)
    plt.ylabel('y',fontsize = 14)
    plt.legend()