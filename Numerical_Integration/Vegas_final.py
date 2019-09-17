import numpy as np
import numpy.random as rn
import matplotlib.pyplot as plt
import time as tm
import sys as sy
import functions as fun
coat_ans  = 0.497661 # answer from the omnipitent machine

#%%
class sbin():
    # set class parameters that can be called from anywhere
    f, rang, r_mag, bk = [None] * 4 # these will be set when vegas is called  
    no_s, Q_tot, Q_sq = [0] * 3 # set counters
    bins, I_store, er_store, x_store = [], [], [], [] # store lists
    f_store = np.zeros(100) # histogram of x_store added to sequentially
    
    def __init__(self,x_1,x_2,v):
        self.min, self.max = x_1, x_2
        self.r = abs(x_1-x_2)
        self.val = v
        self.sum = 0
        self.sq_sum = 0
        
    def sample(self):
        'Give a sample from within this bin'
        x = rn.random() * self.r + self.min # random value between x_1 and x_2 (uniform)
        q = sbin.f(x)/self.val
        self.sum += q
        self.sq_sum += q**2
        sbin.x_store.append(x)

    def info(self):
        'provides the current information of a bin, mainly for debugging'
        print('range ', self.min,':', self.max, ' Value ', self.val)

#%%
def iterate(n, iter_number = 5):
    """Takes n samples from the current pdf then uses them to
    update the pdf. iteration is the current iteration number,
    used to decide if burn out has finished"""
    
    sbin.no_s = n # store the number of iterations
    for i in range(n):
        s = int(rn.random() * len(sbin.bins)) # pick a random bin
        sbin.bins[s].sample() # sample that bin
        
    for i in sbin.bins: # find the totals for all bins
        sbin.Q_tot += i.sum
        sbin.Q_sq += i.sq_sum
        
    er = np.sqrt( (sbin.Q_sq/n - (sbin.Q_tot/n)**2) / (n-1) ) # find error
    update_pdf() # update pdf

    # store values after burn in
    if iter_number >= 5:
        sbin.er_store.append(er)
        sbin.I_store.append(sbin.Q_tot/sbin.no_s)
        
        # store x values sampled in a histogram as you go
        to_hist = np.array(sbin.x_store) * 100/sbin.r_mag
        bc = np.bincount( to_hist.astype(int) )
        sbin.f_store += bc
    
    # reset the totals and x store for the next run
    sbin.Q_tot, sbin.Q_sq  = [0] * 2
    sbin.x_store = []
    
#%%
def update_pdf():    
    """Divides the pdf bins to a fixed total number of sub bins (sbin.bk * 100).
    A bins divison is porportional to how much of the total integral it has.
    Bins are then regrouped so that there is a total of sbin.bk bins again"""
    edges = []
    remainder = 0 # moves fractional bins over to the adjasent bin
    for i in sbin.bins:
        m = ( (sbin.bk * 100) * i.sum / sbin.Q_tot) + remainder # divide a total of 100 * sbin.bk times
        mr = int(np.around(m,0)) # number of times to divide
        remainder = m - mr # remainder for next bin, can be negative
        sub_width = i.r/(mr+1) # divided bins width
        for n in range( mr ):
            edges.append(i.min + (n)*sub_width) # store the start val for new each bin
    
    # if there is a number of bins that can not be evenly group raise
    # if len(edges) % sbin.bk !=0: raise("Error non-divisible bins")
    # should never happen outside debugging so has been removed to improve run time
    
    for i in range( len(sbin.bins) ): # change the bin edges to the new ones
        sbin.bins[i].min = edges[i*100] # group every 100 bins
        try: sbin.bins[i].max = edges[ (i+1)*100 ]
        except IndexError: None # if at final i keep the final bin upper limit
        sbin.bins[i].r = abs(sbin.bins[i].min-sbin.bins[i].max) # update the range
        sbin.bins[i].sum = 0 # reset the sum
    
    for i in sbin.bins: # renormalise all bins, so pdf is uniform in bins
        i.val = 1/(i.r * len(sbin.bins))
   
#%%
def vegas(e, f, lims, total_bins = 200):
    """Implement adaptive importance sampling method Vegas, a uniform pdf of
    bins is edited to be more and more like an optimal pdf for any funciton.
    e is the relative accuracy, f is the integrand, lims are integration range
    [a,b] with a<b, total_bins is the resolution of the generated pdf.
    (note this edits the sbin class vairables and creates total_bins instances)
    returns an array of each iteration and its error, final entry is the
    solution"""
    
    # set parameters in sbins class, allowing them to be called from anywhere
    sbin.f, sbin.rang, sbin.bk = f, lims, total_bins
    sbin.r_mag = abs(lims[0] - lims[1])
    
    #set up the initial uniform bins
    index = np.linspace(lims[0], lims[1], sbin.bk + 1)
    sbin.bins = [ sbin( index[i], index[i+1], 1/sbin.r_mag) 
                  for i in range(sbin.bk) ]
    rn.seed(0) # set the random number seed
    n, c, std_er = 1,1,1
    
    while True: # 10 used as vairence is often low for a small number of points
        n += 1 # increase the iteration counter
        iterate(total_bins * 100, n) # run Vegas once, around 100 samples per bin
        
        if n >= 5:
            # combine the integral estimations weighted by their errors 
            I_i = np.array(sbin.I_store)
            er_i = np.array(sbin.er_store)
            
            #weighted average method
            w = 1/er_i
            I = np.average(I_i, weights=w)
            er = np.average( (I_i - I)**2, weights=w ) 
            std_er = np.sqrt(er/n)   
        
            try: converge = np.vstack((converge, np.array([[I,std_er]])))
            except NameError: converge = np.array([[I,std_er]])
            if n>10 and c > 20: break
            elif std_er < e: c += 1
            
            # display current iteration info
            sy.stdout.write('\rIteration: %.4g     Error %.4g' % (n, std_er) )
                           
    print('\nIntegral is:',I, '+-', std_er, '(epsilon = %f)' % (std_er/e) )    

    return converge # has answer and error as last entery

#%%
def plot_veg(converge, comp_val, fig = 'vegas'):
    """Plot the original funciton and found pdf between the limits
    lim = [a,b]. Histogram is setup to have 100 bins in Vegas.
    Also plots the convergence of the integration over time"""
    
    # plot fuction, pdf and histogram (as a line)
    # x values to plot over, spaced so that histogram works on range given
    x_p = np.arange(sbin.rang[0], sbin.rang[1], sbin.r_mag/100)
    # normalise the histogram plot
    sbin.f_store /= (sbin.f_store.sum() * sbin.r_mag/100)
    
    plt.close(fig), plt.figure(fig)
    plt.plot(x_p, [sbin.f(x) for x in x_p], 'b-', lw = 2, label = 'Funciton')
    plt.plot(x_p, [sbin.f(x)/coat_ans  for x in x_p], 'r-', lw = 2, label = 'Optimal PDF')
    plt.plot(x_p, sbin.f_store ,'k-', label = 'samples histogram')
    plt.legend(), plt.grid()
    plt.title('Optimal PDF & Sampled values'), plt.xlabel('Z'), plt.ylabel('f(Z)')

    # plot convergence
    plt.close(fig + 'converge'), plt.figure(fig + 'converge')
    x_p = converge.shape[0] # x values to plot over
    plt.errorbar( range(x_p), converge[:,0], converge[:,1], fmt = 'none',
                 label = 'Convergence', ecolor='b', capsize=5)
    plt.plot(range(x_p), converge[:,0], 'bo')
    plt.plot( range(x_p), np.ones(x_p)*coat_ans , 'r--', lw = 0.5, label = 'act')
    plt.legend()#, plt.grid()
    plt.title('Covergence of Vegas Algorithm')
    plt.ylabel('Calculated Integral'), plt.xlabel('Iteration')

#%%
def solve_veg(e = 1e-3, f=fun.f_given, lims=[0,2], total_bins = 200):
    """
    Solve vegas for (by defult) the given wave function, plot and validate the
    errors. Call with f, lims and e to do the same for other
    integrads. total_bins is the number of bins in the changing pdf"""
    print('Running...')
    t = tm.clock()
    conv = vegas(e, f, lims, total_bins)
    if abs(coat_ans - conv[-1,0]) < conv[-1,1]:
        print('Agrees with Simpsons Rule Result',
              '(%.4g epsilon)' % (abs(coat_ans - conv[-1,0])/e))
    else:  print('Does not agree with Simpsons Rule Result')
    print('Run time %.4fs' % (tm.clock() - t))
    print( 'Difference from Newton Coats anwer %.4g' % abs(conv[-1,0]-coat_ans),
           '( %.4g epsilon)' % (abs(conv[-1,0]-coat_ans)/e) )
    plot_veg(conv, comp_val = coat_ans)
    return
    
#%%
def val_veg(plot = False):
    """
    Plot = True if convergence plots are wanted.
    Run Vegas integration for test functions specified in functions.py.
    The functions given are: t1 = cos(z), t2 = pi x^3, t3 is a guassian
    with mean 1 and sd 1"""
    # for every function to be tested from functions.py
    for f,lims,e,I in zip(fun.test_funcs, fun.test_ranges,
                          fun.test_accuracies_vegas,
                          fun.test_integrals):
        print('Running ', f.__name__, '...')
        t = tm.clock()
        converge = vegas(e, f, lims, total_bins=100)
        if abs(converge[-1,0]-I) < converge[-1,1]: print('Agrees with true value')
        else: print('Does not agree with true value')
        print( 'True Value is',I,', difference is %.4g' % abs(converge[-1,0]-I),
               '( %.4g epsilon)' % abs((converge[-1,0]-I)/e) )
        print('Run time %.4fs \n' % (tm.clock() - t))
        
        if plot == True: plot_veg(converge, fig = f.__name__, comp_val = I)
    return
    