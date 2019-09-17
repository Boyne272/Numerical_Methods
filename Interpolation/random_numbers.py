import numpy as np
import numpy.random as rn
import matplotlib.pyplot as plt
import time as tm

def reject(C, n):
    'C is Comparison funciton, n is the number to be rejected or not'
    r_num = rn.random() * C(n) # random number between zero and comparison funct at n
    if r_num > pdf(n): # if random number is above pdf at this value   
        return True # reject it
    else:
        return False # keep it

def trans(x):
    """Transformation function (Inverse cumulative distribution function)
    found analytically for 0.65*sin(x)"""
    return np.arccos(1-2*x)

def compar_f(x):
    """Comparison Function for rejcetion method,
    must always greater than or equal to the pdf"""
    return (2/np.pi)*np.sin(x)

def pdf(x):
    'PDF we want the random numbers distributed over'
    return (2/np.pi) * (np.sin(x)**2)


def verify_rn(u_count, u_x, t_count, t_x, r_count, r_x):
    """Check if the comparison funciton is greater than the pdf over the whole range,
    and if the numbers outputted are close to the pdfs."""
    
    if (compar_f(f_range) >= pdf(f_range)).sum() == f_range.shape[0]:
        print (':) Comparison funciton is always greater than pdf')
        
    # needs to be compacted
    t_x = t_x[:-1] # remove the exatra edge from the end
    t_x += (t_x[1]-t_x[0])/2 # center each bin x value
    r_x = r_x[:-1] # remove the exatra edge from the end
    r_x += (r_x[1]-r_x[0])/2 # center each bin x value
    u_x = r_x[:-1] # remove the exatra edge from the end
    u_x += (r_x[1]-r_x[0])/2 # center each bin x value
    
    u_dif = (u_count - n/res)**2
    t_dif = (t_count - 0.5*np.sin(t_x))**2
    r_dif = (r_count - (2/np.pi) * np.sin(r_x)**2)**2
    print('\nThe mean squared difference between distribution generated and expcted are:')
    print('uniform %.5f' % (u_dif.mean()/n), '\n0.5Sin(x) %.5f' % t_dif.mean(), '\n2/pi*sin^2(x) %.5f' % r_dif.mean() )
   
    
def r_method(N):    
    dist = [] # list to store all accepted values
    n_rej = 0 # counter for the number of numbers rejected
    
    t = tm.clock()
    while len(dist) < N: # untill there are N numbers
        n = trans(rn.random()) # generate random number in the range given
        if reject(compar_f, n) == False: # if you dont reject
            dist.append(n)
        else:
            n_rej += 1
    print('Total number rejected:', n_rej)
    print('Time taken to generate', N, 'numbers by rejection is %.4f' % (tm.clock() - t),'s' )
    return dist


rn.seed(10) # set the seed    
n = 5000 # number of random numbers to be generated
res = 150 # used for the resolution of plots and histograms
f_range = np.linspace(0,np.pi,res) # range of the numbers to be generated

uniform = rn.sample(n) # uniform distribution
transform = trans(uniform) # distribution over 0.5*sin(x)
rejection = r_method(n) # distribution over (2/pi)*sin(x)**2

# Plot results
fig, ax = plt.subplots(3,1, figsize = (10,10))
plt.subplots_adjust(hspace = 0.35)

# plot the uniform distriubtion
u = ax[0].hist(uniform, bins = res, color = 'b', lw = .3, ec = 'k', label = 'Distribution') # not normalised
ax[0].plot(np.arange(0,1,1/res), np.ones(res)*n/res, 'r', lw = 2, label = 'Average') # plot epected average
ax[0].set_title('Uniform Distribution')
ax[0].set(xlabel = 'x')
ax[0].legend()

# plot the transormation method
t = ax[1].hist(transform, bins = res, color = 'b', lw = .3, ec = 'k', normed = True, label = 'Distribution')
ax[1].plot(f_range, 0.5*np.sin(f_range), 'r', lw = 2, label = r'$\frac{1}{2}Sin(x)$')
ax[1].set_title('0.5 Sin(x) Distribution - Transformation Method')
ax[1].set(xlabel = 'x')
ax[1].legend()

# plot the rejection method
r = ax[2].hist(rejection, bins = res, color = 'b', lw = .3, ec = 'k', normed = True, label = 'Distribution')
ax[2].plot(f_range, pdf(f_range), 'r', lw = 2, label = r'$\frac{2}{\pi}Sin^{2}(x)$')
ax[2].plot(f_range, compar_f(f_range), 'g--', lw = 2, label = r'$\frac{2}{\pi}Sin(x)$')
ax[2].set_title('n/2 Sin^2(x) Distribution - Rejection Method')
ax[2].set(xlabel = 'x')
ax[2].legend()

verify_rn(u[0],u[1],t[0],t[1],r[0],r[1])