import numpy as np
import copy as cp

def round_error(n):
    """
    n is a number of vairable type you wish to find the accuracy of
    (value of n = 1 for machine accuracy). Returns the smallest value that
    can not be added and stored in the same type.
    """
    n_type = type(n)
    delta = cp.copy(n) # creates a copy of the same type
    
    while n + delta != n: # loops decreasing delta till result rounds down to n
        smallest = cp.copy(delta) # stores the smallest value that can be added
        delta = n_type(delta / 2) # half delta whilst keeping its type
        
    return smallest 

def solve_f():
    """
    Find the machine accuracy for a range of numerical data types all for
    the number one. (n.b. folat-128 not included as it does not work/exist
    on uni machines)
    """
    types = [np.int32, np.int64, np.float16, np.float32, np.float64]
    for i in types:
        print(i, 'has machine accuracy ', round_error( i(1) ) )