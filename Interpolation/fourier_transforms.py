import numpy as np
import numpy.fft as ft
import matplotlib.pyplot as plt
plt.close('all')

def h(t):
    'Signal function'
    if t >= 2 and t <= 4:
        return 5.
    else:
        return 0.

def r(t):
    'Response function'
    return ( 1/np.sqrt(2*np.pi) ) * np.exp(- t**2 /2.)
    return h(t)

# Constants & ranges
N = 2**10 # number of sample points, should be a power of 2
rang = 6
t_range = np.linspace(-rang,rang, N) # time range
dt = t_range[1] - t_range[0]
dw = 2 * np.pi / (N*dt)
f_range = [p * dw for p in np.arange(-N/2, N/2, 1)] # frequency range

h_sample = [h(x) for x in t_range] # sample of each function for reference
r_sample = [r(x) for x in t_range]

r_shift = r_sample[N//2:] + r_sample[0:N//2] # shift the response function to get correct fft 

h_ft = ft.fft( h_sample ) # fft of each funciton
r_ft = ft.fft( r_shift )

f_prod = h_ft * r_ft
t_convol = ft.ifft(f_prod) *dt # should be the convolution in time that is normalised

h_ret = ft.ifft(h_ft) # inverse fft, should be the orginal funcitons
r_ret = ft.ifft(r_ft)

conv = np.convolve(h_sample, r_sample, mode = 'same') *dt # reference convolution

def pl_ft():
    fig_1, ax_1 = plt.subplots(2,2)
    fig_1.suptitle('Functions and ffts')
    
    ax_1[0][0].set_title('Original Functions')
    ax_1[0][0].plot(t_range, r_shift, 'b-', lw = 2)
    ax_1[0][0].set(xlabel = 'time', ylabel = 'r(t)')
    ax_1[1][0].plot(t_range, h_sample, 'b-', lw = 2) 
    ax_1[1][0].set(xlabel = 'time', ylabel = 'h(t)')
    
    ax_1[0][1].set_title('Fourier Transform')
    ax_1[0][1].plot(f_range, np.real(r_ft), 'r-', lw = 1, label = 'Real')
    ax_1[0][1].plot(f_range, np.imag(r_ft), 'b--', lw = 1, label = 'Imag')
    ax_1[0][1].set(xlabel = 'frequency', ylabel = 'F[r(t)]')
    ax_1[0][1].legend()
    ax_1[1][1].plot(f_range, np.real(h_ft), 'r-', lw = 1, label = 'Real')
    ax_1[1][1].plot(f_range, np.imag(h_ft), 'b--', lw = 1, label = 'Imag')
    ax_1[1][1].set(xlabel = 'frequency', ylabel = 'F[h(t)]')
    ax_1[1][1].legend()

def pl_con():
    fig_2, ax_2 = plt.subplots(2,1)
    #fig_2.suptitle('Convolutions')
    
    ax_2[0].set_title('Convolutions') # should be Convolution
    ax_2[0].plot(t_range, np.real(t_convol), 'r-', lw = 1, label = 'Real')
    ax_2[0].plot(t_range, np.imag(t_convol), 'b--', lw = 1, label = 'Imag')
    ax_2[0].set(ylabel = 'Inverse fft product')
    ax_2[0].legend()
    
    #ax_2[1].set_title('Numpy Convolution')
    ax_2[1].plot(t_range, np.real(conv), 'r-', lw = 1, label = 'Real')
    ax_2[1].plot(t_range, np.imag(conv), 'b--', lw = 1, label = 'Imag')
    ax_2[1].set(xlabel = 'time', ylabel = 'Numpy Convolution')
    ax_2[1].legend()

def pl_other():
    fig_3, ax_3 = plt.subplots(2,1, figsize = (10,10))
    fig_3.suptitle('Reference plots')
    
    ax_3[0].set_title('Fourier transform product')
    ax_3[0].plot(f_range, np.real(f_prod), 'r-', lw = 1, label = 'Real')
    ax_3[0].plot(f_range, np.imag(f_prod), 'b--', lw = 1, label = 'Imag')
    ax_3[0].set(xlabel = 'frequency', ylabel = 'F^-1[ r(w) * h(w) ]')
    ax_3[0].legend()
    
    ax_3[1].set_title('Inverse fourier transforms') # should be original functions
    ax_3[1].plot(t_range, np.real(h_ret), 'b-', lw = 1, label = 'h(t)')
    ax_3[1].plot(t_range, np.real(r_ret), 'r-', lw = 1, label = 'r(t)')
    ax_3[1].set(xlabel = 'time', ylabel = 'h(t) and r(t)')
    ax_3[1].legend()
    
pl_con()
pl_ft(),  pl_other()