import numpy as np
import numpy.random as rn
import time as tm

A = np.array([ [2 ,1 ,0 ,0 ,0],
               [3 ,8 ,4 ,0 ,0],
               [0 ,9,20,10 ,0],
               [0, 0,22,51,-25],
               [0, 0, 0,-55,60], ])

b = np.array([ [2],[5],[-4],[8],[9] ] )

class LU_decomp:
    ''
    def __init__(self, Mat):
        """
        Assumed Mat is a square matrix which does not requier Pivoting
        """
        self.Mat = Mat
        self.i, self.j = Mat.shape[0], Mat.shape[1] # keep the dimensions for indexing
        self.Low , self.Upp = self.decomp() # find the upper and lower matrices
        
        # optional commands
        self.f_det() # find Mat determinant
        self.inverse = self.FB_sub( np.identity( self.i ) ) # find Mat inverse
        self.LU_combined = self.Low + self.Upp - np.identity(self.i) # returns combination of L and U as requiered for task
    
    def decomp(self):  
        """
        Find the Upper and Lower Matricies for the Matrix
        Can be slightly optimised by removing the for loop for an array slice
        ******PIVOT ME********
        ***OPT ME*** wont optimise till need to :)
        """      
        L = np.identity(  self.i ) #create an array with diagonal ones
        U = np.zeros([ self.i, self.j ])
        
        for n in range( self.i ): #for every row/column pair
        
            # determine the nth row in U
            for j in np.arange( n, self.i ): # for every column after the diagonal
                U[n,j] = self.Mat[n,j] - np.array( [ L[n,k] * U[k,j] for k in np.arange(0, n)] ).sum()
                
            #determine the nth column in L
            for i in np.arange( n, self.i ):# for all the rows after the diagonal 
                L[i,n] = (1/U[n,n]) * ( self.Mat[i,n] - np.array( [ L[i,k] * U[k,n] for k in np.arange(0, n) ] ).sum() )
        return [L,U]
    
    def f_det(self): #find the determinant
        self.det = 1.
        for i in np.arange( self.i ): # for every diagonal in Upper
            self.det *= self.Upp[i,i]

    def FB_sub(self,b):
        """
        ***OPT ME*******COMMENT ME*****
        """
        x = np.zeros(b.shape)
        for j in np.arange( b.shape[1] ):
            # for every column in b do fb subsitution on b_j column
            y = np.zeros( b.shape[0] )
            x_j = y.copy() # matricies for y and x 
            
            for i in np.arange( y.shape[0] ): # solve Ly = b for y by forward substitution
                y[i] = (1 / self.Low[i,i]) * (b[i,j] - np.array( self.Low[i,:i] * y[:i] ).sum())
                
            for i in np.arange(x_j.shape[0])[::-1]:
                # solve Ux = y for x by backward substitution, being sure to solve for i = N-1 first
                s = self.Upp[ i, i+1:x_j.shape[0] ] *x_j[ i+1:x_j.shape[0] ]
                x_j[i] = (1/self.Upp[i,i]) * ( y[i] - np.array( s ).sum() )
                
            x[:,j] = x_j # take this column solution and put it in the total matrix
        return x
    
    def check(self):
        """
        Tests if the LU product is original matrix within float errors and
        if the Mat determinant is equal to the product of the Upper matrix diagonals 
        and if  inv_Mat * Mat = Identity  
        """
        LU = np.matmul( self.Low, self.Upp)
        C = abs(self.Mat-LU) # Mat - LU elements should be zero within rounding errors
        if C.sum() < 1E-5: # arbitary 1E-5 used as this is the sum of rounding errors
            print('Correct :)   L*U = Mat')
        else:
            print('Error :(   L*U != Mat within 1E-5 total rounding error')
         
        det = np.linalg.det(self.Mat) # find the determinant with numpy
        
        #print('LU_decomp det = ', self.det, '\nNumpy det = ', det)
        
        if abs( self.det - det ) <= abs(det/1E5):
        # a fravtion of det is used as rounding errors are relative to size of number
            print('Correct :)   determinant is right')
        else: 
            print('Error :(   determinant is Wrong or too big (check if det A = inf)')
            
        I = np.matmul( self.Mat , self.inverse)
        if (I.sum() - I.shape[0]) < 1E-5:
            print('Correct :)   mat * mat^-1 = I ')
        else:
            print('Error :(   Mat * Mat^-1 != I within 1E-5 total rounding error')
        print('All within 1E-5 sum of rounding errors (or 1E-5 fraction for det)')

def validate():
    ''
    # test arrays, a random array is used as defult
    #T = np.array( [ [4.,3.],[6.,3.] ] )
    #T = np.array([ [1.,2.,4.] , [3.,8.,14.] , [2.,6.,13.] ])
    T = rn.random_sample([200,200]) * 100
    t = tm.clock()
    Test = LU_decomp(T) 
    print('\nTime to solve test: %.5fs' % ( tm.clock() - t ) )
    Test.check()
    return Test

Sol = LU_decomp(A)
x = Sol.FB_sub(b)
print('\nsolved Ax=b with A:\n',A,'\nb:\n',b,'\nsolution x:\n',x)
#check
Sol.check() 