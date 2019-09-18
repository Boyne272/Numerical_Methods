# ACSE-5  Project  3  -  Positive  Definite  Matrix  Solver

Richard Boyne (team mr fantabulous)
	- Notes about the file structure
	- How to use for own matrices
	- How the generate matrix method works (sudo code???)
	- User interface for main.cpp (basically a bunch of examples)
  
## Code Layout
There are four main secions of code:
- matrix.h, matrix.cpp hold all the methods to make a matrix and solve them
- generate_matrix.cpp holds the methods for creating a random diagonally dominant matrix with a given sparcity
- menu.h is from https://github.com/Boyne272/menu_cpp (create by me) to give a user interface
- main.cpp has a selection of coding examples of how one should use the matrix classes. It also includes some functions for timing the code and saving it which was used to generate the report. All of this is wrapped in the menu user interface so simply running the code allows one to see these examples in action.

## Solving your own matrix
First one needs to initialise the matrix by calling the appropirate constructor:
- dense with no initial values, just give the dimensions
- dense with all elements set to an intitial value
- sparse with dimensions and number of non-zero entries to be entered

Second values need to be set by calling the appropirate set method:
- dense matrix this can be done by calling set with one value at a time and this will cycle through the matrix and set each value
- dense matrix can also be called with a column, row and value (use the 2nd constructor to ensure all values are intialised)
- sparse matrix must be set in the correct order with elements as they occure going along each row

Third the vector b must be set the same as above but with the column dimension set to 1

Finally call the gauss-seidel method of the matrix with the vector b and a number of iterations to solve for

## Random matrix generation
This is done by first choosing an appopirate number of off-diagonal elements to ensure the desired sparsity. These values are generated with the rand() function, which is seeded first with srand(), then modulused by a user given value to ensure it remains small. Once generated they are placed in the non-diagonal matrix enteries at regular intervals to span the whole matrix. Finally the diagonals are calculated by summing all the elements in the row. 

This algorithm is implemeneted in both rand_PDM_dense and rand_PDM_sparse functions with a random seed set such that the matrices produced are identical and can be compared. 
