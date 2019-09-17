import numpy as np
import numpy.random as rn
import matplotlib.pyplot as plt
import time as tm
import sys as sy
import functions as fun
coat_ans = 0.497661 # agreed answer from the newton coats method 

#%%
def trans_int(f, pdf, trans, lims, e, consec = 50, sd = 0):
    """
    Monte Carlo intergation using importance smapling from a given pdf, sampled
    with tranformation method.
    
    f is the integrand (must be callable), pdf is the function to distribute
    samples on normalised on the integration range, trans is assoicated CDF.
    lim is [a,b] to integrate over (a<b), consec determines how strict the
    convergence criteria is, sd is the seed for the random number generator
    
    returns array of integral and error values for each iteration (final 
    entery is solution) and a histogram of all samples taken in an array"""
    rn.seed(sd) #seed the generator
    n, c = 1,0  #counter of samples taken & break counter
    Q_sum, Q_sq_sum = 0,0
    x_samples = [] # temp store of x samples
    hist_samples = np.zeros(100)  # histogram of x_samples added to sequentially
    hist_width = 100 / abs(lims[1] - lims[0]) 
    
    while(True):
        new_x = trans(rn.rand()) # random number from distributed over pdf
        x_samples.append(new_x)
        Q = f(new_x)/pdf(new_x) # find new Q value
        Q_sum +=  Q # add to the total
        Q_sq_sum += Q**2
        n += 1 #increase the counter
        
        # every 1000 steps after the burn inn
        if n % 1000 == 0 and n>5000:
            I = Q_sum / n # new estimate guess
            
            # histogram all samples to reduce storage space
            to_hist = np.array(x_samples) * hist_width # scale to a list from 0 to 100
            bc = np.bincount( to_hist.astype(int) ) # count the enteries rounded down
            hist_samples += bc # add these counts to store
            
            #find the standard error from vairience
            error = np.sqrt(( Q_sq_sum/n -  (Q_sum/n)**2) /n)
            try:
                converge = np.vstack(( converge,[[I,error]] )) # store the values
                rel = abs( 1. - converge[-1,0] / converge[-2,0]) # relative accuracy
                if c >= consec: break # if accurate enough consec times
                elif rel < e: c += 1
                #if error < e*I: break
            except NameError: converge = np.array([[I,error]]) # create store first time
            
            sy.stdout.write('\rIteration: %.4g     Error %.4g' % (n, error) )
            
    print('\nTransformation Monte Carlo found Integral = %.6g' % I,
          '+- %.6g' % error, 'after', n, 'samples')
    
    return converge, hist_samples

#%%
def plot_tran(sam, con, f, pdf, lims, comp_val, fig = 'tran'):
    """Plot the funciton, pdf and normalized histogram 
    of samples used in the algorithmbetween the limits lim = [a,b].
    Also plot the convergence of the algorithm over iterations.
    sam is the histogramed samples, con is an array of I values
    after each iteration with errors (generated from trans_int), comp_val
    is line plotted on convergence graph if given"""
    # plot the funciton, pdf and sampled values
    plt.close(fig), plt.figure(fig)
    x_p = np.arange(lims[0], lims[1], abs(lims[0]-lims[1])/100 ) # x values to plot over
    plt.plot(x_p, [f(x) for x in x_p], 'b-', lw = 2, label = 'Funciton')
    plt.plot(x_p, [pdf(x) for x in x_p], 'r-', lw = 2, label = 'PDF')
    # plot a normalised histogram of samples if given
    sam *= 100/( sam.sum() * abs(lims[0]-lims[1]) ) # normalise counts over range
    plt.plot(x_p, sam, 'k', label = 'Samples')
    plt.title('Integrand and Sampled Points Histogram')
    plt.xlabel('Z'), plt.ylabel('f(Z)')
    plt.legend(), plt.grid()
    
    # plot convergence
    plt.close(fig+' con'), plt.figure(fig+' con')
    x_p = np.arange(con.shape[0])[::2]
    plt.errorbar(x_p, con[:,0][::2], yerr= con[:,1][::2], label = 'Est. Error',
                 capsize = 5, ecolor='r', elinewidth=1, fmt='none')
    plt.plot(x_p, con[:,0][::2], 'ko', label = 'Convergence')
    plt.plot(x_p, np.ones(con.shape[0])[::2]*comp_val, 'b--', lw = 1.5, label = 'Comp_Val')
    plt.title('Covergence of Monte Carlo Transformation')
    plt.xlabel('iteration'), plt.ylabel('Calculated Integral')
    plt.grid(), plt.legend()
    
#%%
def solve_tran(f = fun.f_given, pdf = fun.pdf_given, trans = fun.trans_given,
               lims = np.array([0.,2.]), e = 1e-4):
    """
    Solve tranformation integration for (by defult) the given wave function,
    using the given pdf and its associated CDF plot. Call with f, pdf, trans,
    lims and e to do the same for other integrads"""
    print('Running...')
    t = tm.clock()
    
    converge, samples = trans_int(f, pdf, trans, lims, e)
    if abs(converge[-1,0]-coat_ans) < converge[-1,1]: # if diff < error
        print('Agrees with Simpsons Rule Result')
    else: print('Does not agree with simpsons result')
    
    print('Run time %.4fs' % (tm.clock() - t))
    plot_tran(samples, converge, f, pdf, lims, coat_ans)
    return
     
#%%
def val_tran(plot = False): 
    """
    Plot = True if convergence plots are wanted.
    Run Transformation integration for test functions specified in functions.py
    using the pdfs and cdfs also given there.
    The functions given are: t1 = cos(z), t2 = pi x^3, t3 is a guassian
    with mean 1 and sd 1"""
    # for every function to be tested from functions.py
    for f,pdf,trans,lims,e,I in zip(fun.test_funcs, fun.test_pdfs,  
                                   fun.test_trans, fun.test_ranges,
                                   fun.test_accuracies, fun.test_integrals):
        print('Running ', f.__name__, '...')
        t = tm.clock()
        converge, samples = trans_int(f, pdf, trans, lims, e)
        if abs(converge[-1,0]-I) < converge[-1,1]: print('Agrees with true value')
        else: print('Does not agree with true value')
        print( 'True Value is',I,', difference is %.4g' % abs(converge[-1,0]-I),
               '( %.4g epsilon)' % abs((converge[-1,0]-I)/e) )
        
        print('Run time %.4fs \n' % (tm.clock() - t))
        if plot == True: plot_tran(samples, converge, f, pdf, lims,
                                   comp_val = I, fig = f.__name__)
    return

#%%
def plot_times():
    'plot the iterations taken with relative accuracy (data collected seperatley)'
    plt.figure('Iterations')
    plt.title('Requiered Iterations for Monte Carlo')
    x = [3,4,5,6]
    y_uni = [.74000, 1.08000, 3.24000, 10.05000]
    y_pdf = [.60000, 1.86000, 5.98000, 17.76000]
    plt.plot(x,y_uni,'bo', label = 'Uniform PDF')
    plt.plot(x,y_uni,'b--')
    plt.plot(x,y_pdf,'ro', label = 'Linear PDF')
    plt.plot(x,y_pdf,'r--')
    plt.xlabel(r'Relative Accuracy Exponent Magnitude')
    plt.ylabel(r'Iterations Requiered ($\times 10^{5}$)')
    plt.legend()
    return