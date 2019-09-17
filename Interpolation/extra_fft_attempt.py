import numpy as np
import matplotlib.pyplot as plt

# constants
dig = 5 # the number of samples power
N = 2**dig # number of points =
r = 2 # pm range
dt = 2 * r / N
dw = 2 * np.pi / (N*dt)

# ranges
pr = np.arange(-N/2, N/2, 1)
f = pr * dw
x = np.arange(-r, r , dt)

#function
def h(t):
    'Signal function'
    if t >= -0.5 and t <= 0.5:
        return 5.
    else:
        return 0.
    #return np.sin(t * np.pi/2)

# define dft
def discrete_ft(s, p):
    sum = 0 + 0j
    for n in range( s.shape[0] ):
        sum += s[n] * np.exp( 1j * 2 * np.pi * p * n / N )
    return sum

def fft(f):
    fN = f.shape[0]
    if fN == 1:
        return np.array([f[0]])
    f_even = fft(f[::2])
    f_odd = fft(f[1::2])
    s = [np.exp(1j * 2 * np.pi * p/fN) for p in range(fN//2)]
    p_1 = [f_even[p] + s[p] * f_odd[p] for p in range(fN//2)]
    p_2 = [f_even[p] - s[p] * f_odd[p] for p in range(fN//2)]
    #if fN == N: # fixes frequency shift
    #    return np.array(p_2+p_1)
    return np.array(p_1 + p_2)

def reverse(n,d):
    #TESTED WORKS
    'n is the number and d is the digits you want to reverse over'
    bin_n = bin(n) # convert number to binary
    while len(bin_n) < d+2: # pad out smaller nymbers with zeros 
        bin_n = bin_n[:2] + '0' + bin_n[2:] #so they make right bigger numbers
    a = bin_n[0:2] #include the first two characters that are flags
    for i in bin_n[:1:-1]: # make a reverse number in string
        a += i
    return int(a,base=0) # return string to binary    

# function sample
sample = np.array([h(i) for i in x])
# find my dft
dft = [discrete_ft(sample, q) for q in pr]

# my fft
out = fft(sample)
print('done')
mft = np.zeros( N ).astype(complex)
for o in range( N ):
    rev = reverse(o,dig)
    mft[rev] = out[o]    

# Plot original function
fig_2, ax = plt.subplots(2,2, figsize=(12,8))
ax[0][0].plot(x,sample,'ko')
ax[0][0].plot(x,sample,'k-')

# Plot my dft
ax[0][1].plot(f, np.real(dft),'bo')
ax[0][1].plot(f, np.real(dft),'b-')
#ax[1].plot(f, np.imag(dft),'ro')
#ax[1].plot(f, np.imag(dft),'r-')

#solve by numpy fft
ax[1][0].plot( f, np.real(np.fft.fft([h(i) for i in x])), 'bo' )
ax[1][0].plot( f, np.real(np.fft.fft([h(i) for i in x])), 'b-' )
#ax[1][0].plot( f, np.imag(np.fft.fft([h(i) for i in x])), 'ro' )
#ax[1][0].plot( f, np.imag(np.fft.fft([h(i) for i in x])), 'r-' )

# my fft
ax[1][1].plot(f, np.real(out), 'bo')
ax[1][1].plot(f, np.real(out), 'b-')
#ax[1][1].plot(f, np.imag(mft), 'bo')
#ax[1][1].plot(f, np.imag(mft), 'b-')

ax[0][0].set_title('Sample Function')
ax[0][1].set_title('Real/Imaginary dft')
ax[1][0].set_title('Numpy fft')
ax[1][1].set_title('My fft')