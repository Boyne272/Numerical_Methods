import numpy as np
import matplotlib.pyplot as plt
import time as tm
import functions as fun

#%%
def trapezoidal(f, lims, e):
    """
    Implement Trapizodial rule upto to relative accuracy e.
    f is integrand, lims is an array [lower limit, upper limit]
    returns array of integral and error values for each iteration (final 
    entery is solution)"""
          
    if not callable(f) or not all([ isinstance(i,float) for i in (lims[0],lims[1],e) ]):
        raise TypeError('One of the entered parameters is of wrong type,',
                        '(numbers must be floats) see doc string')
    
    h = abs(lims[1]-lims[0]) # initial bin width
    T = 0.5 * h * ( f(lims[0]) + f(lims[1]) ) # Trapezium first guess
    b = 1 # number of bins
    
    while True: # this way round so it alwys runs at least once
        h /= 2 # half the sample width
        b *= 2 # number of bins
        
        new_p = np.arange(lims[0]+h, lims[1], h*2) # every new point
        new_f = np.array( [f(x) for x in new_p] ) # f value at the new points 
        T = 0.5 * T + h * new_f.sum() # find new trapezium value
        
        try: # store the trapezium value, find rel error and see if you stop
            if converge[-1,0] == 0.: raise ZeroDivisionError(
            'iteration found itegral as zero. Check itegral is not symetric')
            rel = abs( 1 - T/converge[-1,0]) 
            converge = np.vstack([ converge, np.array([[T,rel*T]]) ]) # store val and error
            if ( rel < e): break # stop if accurate enough
        except NameError: converge = np.array([[T,0.1*T]]) # # if first time create the store, guess 10% error        
    
    print("Trapezoidal Rule found", T, '+-', rel*T, 'using', b, 'bins')
    return converge

#%%
def simpsons(f, lims, e):
    """
    Implement Simpsons rule upto to relative accuracy e.
    f is integrand, lims is an array [lower limit, upper limit]
    returns array of integral and error values for each iteration (final 
                                                                   entery is solution)"""
    
    if not callable(f) or not all([ isinstance(i,float) for i in (lims[0],lims[1],e) ]):
        raise TypeError('One of the entered parameters is of wrong type,',
                        '(numbers must be floats) see doc string')
    
    h = abs(lims[1]-lims[0]) # initial bin width
    T_1 = 0.5 * h * ( f(lims[0]) + f(lims[1]) ) # Trapezium first guess
    b = 1 # number of bins
    
    while True: # always runs at least once
        h /= 2 # half the bin width
        b *= 2 # double the bin counter
        
        new_p = np.arange(lims[0]+h, lims[1], h*2) # every new point
        new_f = np.array( [f(x) for x in new_p] ) # f value at the new points
        T_0 = T_1 # keep previous Trapezium value 
        T_1 = 0.5 * T_0 + h * new_f.sum() # find new trapezium value
        S_j = (4/3) * T_1 - (1/3) * T_0 # new simpsons value
        
        try: # store the simpson value, find rel error and see if you stop
            if converge[-1,0] == 0.: raise ZeroDivisionError(
            'iteration found itegral as zero. Check itegral is not symetric')
            rel = abs( 1 - S_j/converge[-1,0])
            converge = np.vstack([ converge, np.array([[S_j,rel*S_j]]) ]) # store val and error
            if ( rel < e): break # stop if accurate enough
        except NameError: converge = np.array([[S_j,0.1*S_j]]) # if first time create the store, guess 10% error
    
    print("Simpson's Rule found",S_j, '+-', rel*S_j, 'using', b, 'bins')
    return converge

#%%
def error_check_ts(T_con, S_con, fig = 2):
    """
    Plot a log graph with magnitude of errors for both methods along with
    expected error orders, h^2 and h^4. Since there is an unknown constant
    dependent on the second derivative of the functions used the lines should
    have the same gradient (in log space) but be of set"""
    plt.close(fig), plt.figure(fig)
    
    x_p1 = np.arange(T_con.shape[0]) # x to plot over
    plt.plot(x_p1, T_con[:,1], 'bo', label = 'Trapez Error')
    h_sq = ( 1/2**(x_p1+1) )**2 # h^2, +1 as the 1 bin value not stored
    plt.plot(x_p1, h_sq, 'b--', label = r'$\mathcal{O} h^{2}$')
    
    x_p2 = np.arange(S_con.shape[0]) # different x to plot over
    plt.plot(x_p2, S_con[:,1], 'ro', label = 'Simpson Error')
    h_4 = ( 1/2**(x_p2+1) )**4 # h^4
    plt.plot(x_p2, h_4, 'r--', label = r'$\mathcal{O} h^{4}$')
    
    plt.legend(), plt.grid(), plt.yscale('log')
    plt.title('Error Convergence'), plt.xlabel('Iteration'), plt.ylabel('Error Magnitude')
    return

#%% 
def plot_ts(T_con, S_con, comp_val = None, fig=1):
    """
    Plot the convergence of Trapezodial and Simpsons integration methods,
    fig is the label of the figure to allow for multiple successive calls"""
    
    plt.close(fig), plt.figure(fig)
    #plot trapezium
    x_p1 = np.arange(T_con.shape[0]) # x to plot over
    plt.errorbar(x_p1, T_con[:,0], yerr = T_con[:,1], capsize= 5, ecolor='b',
                 fmt = 'none', label = 'Trapez')
    plt.plot(x_p1, T_con[:,0], 'bo')
    #plot sumpsons
    x_p2 = np.arange(0.2, S_con.shape[0]+0.2)  # x to plot over, ofset by 0.2 to make the plot clear
    plt.errorbar(x_p2, S_con[:,0], yerr = S_con[:,1], capsize= 5, ecolor='r',
                 fmt = 'none', label = 'Simpson')
    plt.plot(x_p2, S_con[:,0], 'ro', label = 'Simpson')
    #plot the line for true value if given 
    x_p3 = np.arange( max(T_con.shape[0], S_con.shape[0]) )  # x to plot over
    if comp_val != None: plt.plot(x_p3, np.ones(x_p3.shape[0])*comp_val,
                                  'g-', label = 'Comparison Value')
        
    plt.legend()#, plt.yscale('log')
    plt.title('Covergence of Simpsons & Trapezium Rules')
    plt.ylabel('Calculated Integral'), plt.xlabel('Iteration')
     
#%%
def solve_ts(f = fun.f_given, lims = np.array([0.,2.]), e = 1e-6):
    """
    Solve Trapezium and Simpsons for (by defult) the given wave function, plot
    and validate the errors. Call with f, lims and e to do the same for other
    integrads"""
    print('Running...')
    t = tm.clock()
    T = trapezoidal(f, lims, e)
    S = simpsons(f, lims, e)
    print('Run time %.4fs' % (tm.clock() - t))
    print( 'Difference between found integrals %.4g' % abs(T[-1,0]-S[-1,0]),
           '( %.4g epsilon)' % abs((T[-1,0]-S[-1,0])/e) )
    error_check_ts(T,S)
    plot_ts(T, S, S[-1,0])
    return

#%%
def val_ts(plot = False, error = False):
    """
    Plot = True if convergence plots are wanted, error = True if error
    validation is wanted.
    Run both simpsons and trapezium for test functions specified in functions.pyand test the zero integral
    estimate is working with an anti-symetric function.
    The defult functions given are: t1 = cos(z), t2 = pi x^3, t3 is a guassian
    with mean 1 and sd 1, symetric function sin(x)"""
    # for every function to be tested from functions.py
    for f,lims,e,I in zip(fun.test_funcs, fun.test_ranges,
                          fun.test_accuracies, fun.test_integrals):
        print('Running ', f.__name__, '...')
        t = tm.clock()
        
        T = trapezoidal(f, lims, e)
        print( 'True Value is',I,', difference is %.4g' % abs(T[-1,0]-I),
               '( %.4g epsilon)' % abs((T[-1,0]-I)/e) )
        
        S = simpsons(f, lims, e)
        print( 'True Value is',I,', difference is %.4g' % abs(S[-1,0]-I),
               '( %.4g epsilon)' % abs((S[-1,0]-I)/e) )
        
        print('Run time %.4fs \n' % (tm.clock() - t))
        if plot == True: plot_ts(T, S, I, fig = f.__name__)
        if error == True: error_check_ts(T, S, fig = f.__name__ + ' error')
        
    # test the zero exception with test_f2 (symetic function)
    try:
        trapezoidal(fun.test_symetric, [-1.,1.], 1e-4) # symetric function on symetric interval
        print('\nError: Trapezium Zero Exception not working\n')
    except ZeroDivisionError: print('Trapezium zero exception works')
    try:
        simpsons(fun.test_symetric, [-1.,1.], 1e-4)
        print('\nError: Simpsons Zero Exception not working\n')
    except ZeroDivisionError: print('Simpsons zero exception works')
    return
