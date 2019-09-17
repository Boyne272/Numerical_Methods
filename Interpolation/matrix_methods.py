import numpy as np
import numpy.random as rn
import time as tm

class LU_decomp:
    """
    Create an instance for each Matrix with methods to find and store
    its original values, upper and lower triangular matrix, its determinant
    and inverse.'
    """
    def __init__(self, Mat):
        """
        Pass in matrix Mat, assumed Mat is a square matrix with no zeros on the
        diagonal (i.e. no pivoting requered) and is not singular.
        """
        self.Mat = Mat
        self.i = Mat.shape[0] # keep the dimension for indexing
        self.Low , self.Upp = self.decomp() # find the upper and lower matrices
        
        # optional commands
        self.f_det() # find determinant
        self.inverse = FB_sub( self.Low , self.Upp, np.identity( self.i ) ) 
        # find the combination of L and U as requiered
        self.LU_combined = self.Low + self.Upp - np.identity(self.i)
        
    def decomp(self):  
        """
        Find the Upper and Lower Matricies for the Matrix.  
        """      
        L = np.identity(  self.i ) # create an array with diagonal ones
        U = np.zeros([ self.i, self.i ])
        
        for n in range( self.i ): # for every row/column pair
        
            # determine the nth row in U
            for j in np.arange( n, self.i ): # for every column after the diagonal
                U[n,j] = self.Mat[n,j] - np.array( [ L[n,k] * U[k,j] for k in np.arange(0, n)] ).sum()
                
            #determine the nth column in L
            for i in np.arange( n, self.i ):# for all the rows after the diagonal 
                L[i,n] = (1/U[n,n]) * ( self.Mat[i,n] - np.array( [ L[i,k] * U[k,n] for k in np.arange(0, n) ] ).sum() )
        return [L,U]                  
    
    def f_det(self):
        """
        find the determinant from the product of upper diagonal matrixcies.
        Be cautious of large matricies with determinants that will overflow
        a float 64.
        """
        self.det = 1.
        for i in np.arange( self.i ): # for every diagonal in Upper
            self.det *= self.Upp[i,i]

    def check(self):
        """
        Tests if the LU product is original matrix within rounding errors and
        if the Mat determinant is equal to the product of the Upper matrix diagonals 
        and if  inv_Mat * Mat = Identity.
        Within rounding errors means the sum of all rounding errors are less
        than 1E-5 (or 0.001% for det)
        """
        print('\nChecking...\n')
        
        LU = np.matmul( self.Low, self.Upp)
        C = abs(self.Mat-LU) # Mat - LU elements should be zero within rounding errors
        if C.sum() < 1E-5: # arbitary 1E-5 used as this is the sum of rounding errors
            print('Correct :)   L*U = Mat')
        else:
            print('Error :(   L*U != Mat within 1E-5 total rounding error')
         
        det = np.linalg.det(self.Mat) # find the determinant with numpy
        print(' LU_decomp det = ', self.det, '\n Numpy det = ', det)
        
        if (self.det / det ) - 1 <= 1E5:
        # a fravtion is used as rounding errors are relative to size of number
            print('Correct :)   determinant is in agreement (check above values are sensible)') # numpy det may also have issues
        else: 
            print('Error :(   determinant is Wrong with rounding error of ', 
                  '0.001% or too big (check if det A = inf)')
            
        I = np.matmul( self.Mat , self.inverse) # should be intentity
        if (I.sum() - I.shape[0]) < 1E-5: # so sum should equal the number of rows
            print('Correct :)   mat * mat^-1 = I ')
        else:
            print('Error :(   Mat * Mat^-1 != I within 1E-5 total rounding error')
        print('\nAll checked within 1E-5 sum of rounding errors (or 1E-5 fraction for det)')

def FB_sub(L, U, b):
    """
    Solve the Matrix equation Mat*x = b, via L*y = b and U*x=y for x. 
    b can be a matrix.
    """
    x = np.zeros(b.shape)
    # for every column in b do fb subsitution on b_j column
    for j in np.arange( b.shape[1] ):
            
        y = np.zeros( b.shape[0] )
        x_j = y.copy() # set up matricies for y and x 
        
        # solve Ly = b for y by forward substitution
        for i in np.arange( y.shape[0] ):
            y[i] = (1 / L[i,i]) * (b[i,j] - np.array( L[i,:i] * y[:i] ).sum())
            
        # solve Ux = y for x by backward substitution, being sure to solve for i = N-1 first
        for i in np.arange(x_j.shape[0])[::-1]:
            # s defimimed then immidiatly used for clarity of code, can combined next rows
            s = U[ i, i+1:x_j.shape[0] ] *x_j[ i+1:x_j.shape[0] ]
            x_j[i] = (1/U[i,i]) * ( y[i] - np.array( s ).sum() )
            
        x[:,j] = x_j # take this column solution and put it in the total matrix
    return x

def validate_LU(size=100, r=10, s=0):
    """
    Runs LU_decomp and forward backward substitution for a random (with seed s) 
    test array of dimension of n,m = size and values ranging from -r to r then
    checks the found LU matrices, determinant and inverse.
    Chances of a singular or pivot requiereing random matrix is negligable.
    Returns the class instance used for the test.
    """
    # other test arrays
    #T = np.array( [ [4.,3.],[6.,3.] ] )
    #T = np.array([ [1.,2.,4.] , [3.,8.,14.] , [2.,6.,13.] ])
    rn.seed(0) # seed the random generator
    T = rn.random_sample([size,size]) * r * rn.choice([-1,1])
    t = tm.clock()
    Test = LU_decomp(T) 
    print('\nTime to solve test %.5fs of size' % ( tm.clock() - t ), size)
    Test.check()
    return Test

def solve_LU():
    """
    Solves the given problems in project A task 2 as requiered and checks
    the solutions.
    """
    
    A = np.array([ [2 ,1 ,0 ,0 ,0], # set up the matricies
                   [3 ,8 ,4 ,0 ,0],
                   [0 ,9,20,10 ,0],
                   [0, 0,22,51,-25],
                   [0, 0, 0,-55,60], ])
    b = np.array([ [2],[5],[-4],[8],[9] ] )
    
    Sol = LU_decomp(A) # solve for LU
    print('\nA:\n',A,'\nDecomposes into LU cobination matrix:\n',
          Sol.LU_combined,'\nand has determinant %.5f' % Sol.det)
    
    x = FB_sub(Sol.Low, Sol.Upp, b) # solve A*x = b
    print('\nsolved A*x=b with, b:\n',b,'\nto find x:\n',x)
    print('\nsolved for inverse of A:\n',Sol.inverse)
    Sol.check() # check solution