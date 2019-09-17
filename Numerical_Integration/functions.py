import numpy as np
import matplotlib.pyplot as plt

#%%
def f_given(z):
    """
    Integrand to be interagted over range [a,b]
    (assumed continuous over this range)""" 
    return (np.abs( np.exp(- z**2 / 2) / np.sqrt(np.sqrt(np.pi)) ))**2

def pdf_given(z):
    """
    pdf for metropolis algorithm sampling, must have simular shape
    to func, normalised and non-zero over the region [0,2]"""
    return -0.48*z + 0.98 # recommended comparison function

def trans_given(z):
    """
    Integrated and inverted of given_pdf for transformation method"""
    a,b = 0.24, 0.98 # y = -ax + b 
    return (b - np.sqrt(b*b- 4*a*z )) / (2*a)

#%%
def uni_pdf(z):
    """
    pdf for metropolis algorithm sampling, assumed the integral limits are 0,2"""
    return 0.5

def trans_uni(z):
    """
    Integrated and inverted of uniform pdf for transformation method"""
    return 2 * z

#%%
def plot_transform(pdf = pdf_given, trans = trans_given):
    """
    Plot given pdf, associated Cumulative Distribution Function
    and histogram of 10,000 random samples over the pdf over the range [0,2]"""
    x_p = np.linspace(0,2,100)
    plt.close(0), plt.figure(0)
    plt.plot(x_p, [pdf(x) for x in x_p], 'r-',label = 'PDF')
    #plt.plot(x_p, [-(0.24 * y**2) + 0.98*y for y in x_p], 'b--',label = 'CI')
    x_p = np.linspace(0,1,100)
    plt.plot(x_p, [trans(x) for x in x_p], 'g-',label = 'CDF')
    plt.hist( [trans(r) for r in np.random.rand(10000)], bins= 100, fc = 'w', ec = 'k',
             normed = True, label = 'Samples')
    plt.legend(), plt.xlabel('x'), plt.ylabel('y')
    plt.title('Transformation Sample Method Plots')
    
#%%
# Validation functions
def test_f1(z):
    """
    first test function, cos(z). The integral on [0,pi/2] is 1"""
    return np.sin(z)

def test_pdf1(z):
    """
    associated pdf for test_f1 over the range [0,pi/2]"""
    return 8*z/(np.pi**2)

def tran_pdf1(z):
    """
    Integrated and inverted of test_pdf1 for transformation method"""
    return np.sqrt(z*(np.pi**2)/4)

def test_f2(z):
    """
    second test function, cubic function times by pi for complexity. integral
    on [0,+2] = 4/pi"""
    return z**3 / np.pi

def test_pdf2(z):
    """
    associated pdf for test_f2 over the range [0,+2]"""
    return 3 * z**2 / 8

def tran_pdf2(z):
    """
    Integrated and inverted of test_pdf2 for transformation method"""
    return 2 * np.cbrt(z)

def test_f3(z):
    """
    third test function, gaussian distribution centered at z = 1
    with standard deviation 1. The integral on [0,2] is 0.6827"""
    return np.exp(-0.5 * (z-1)**2) / np.sqrt(2*np.pi)   

def test_pdf3(z):
    """
    associated pdf for test_f3 over the range [0,2]"""
    if z <=1: return z
    else: return -z + 2

def tran_pdf3(z):
    """
    Integrated and inverted of test_pdf3 for transformation method
    This is done by eperatly integrating each line and then transforming
    from their (integral is always positive so this works)"""
    if z <0.5: return np.sqrt(2*z)
    else: return 2 - np.sqrt(2-2*z) # if 0.5<z<1
    
def test_symetric(z):
    """
    Used to check the zero integral exception is working in valication"""
    return np.sin(z)

# list functions to be tested, the test ranges and expected integrals
test_funcs = [test_f1, test_f2, test_f3 ]
test_pdfs = [test_pdf1, test_pdf2, test_pdf3]
test_trans = [tran_pdf1, tran_pdf2, tran_pdf3]
test_ranges = np.array([ [0,np.pi/2], [0,2], [0,2] ])
test_integrals = [1., 4./np.pi, 0.6827]
test_accuracies = [1e-4, 1e-4, 1e-4]
test_accuracies_vegas = [1e-3, 1e-2, 5e-2] # must be lower for reasonable run time