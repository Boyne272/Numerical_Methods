#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using namespace std;

template <class T>
void printStuffMatrix(Matrix<T>& input_mat)
{
   cout << "Hello!";
}

int main()
{

   int rows = 5;
   int cols = 5;

   // Testing our matrix class
   // CHANGE THIS TO BE TEMPLATED WITH INTS 
   auto *dense_mat = new Matrix<double>(rows, cols, true);
   // Now we need to go and fill our matrix
   // Let's just fill it with random values for now
   for (int i = 0; i < rows * cols; i++)
   {
      dense_mat->values[i] = i;
   }
   // Now let's test printing our matrix
   dense_mat->printValues();
   // Now let's test printing our matrix with our other function
   dense_mat->printMatrix();   

   // Don't forget to call delete (ie the destructor) once we're done, otherwise 
   // we will leak memory from the new called inside our matrix
   delete dense_mat;

   // ~~~~~~~~~~~~~~
   // Now let's try our new sparse class
   // ~~~~~~~~~~~~~~
   int nnzs = 5;
   auto *sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
   // Let's just fill it with values on the diagonals
   for (int i = 0; i < nnzs; i++)
   {
      sparse_mat->values[i] = i;
      sparse_mat->col_index[i] = i;
      sparse_mat->row_position[i] = i;
   }

   sparse_mat->row_position[nnzs] = nnzs;

   // Now let's print it
   sparse_mat->printMatrix();

   // Now let's test our matvec
   double *input = new double[cols];
   double *output = new double[rows];
   for (int i = 0; i < cols; i++)
   {
      input[i] = i;
   }
   // Let's compute our matrix-vector product
   sparse_mat->matVecMult(input, output);
   // Let's print the results
   cout << "Output of our matvec" << endl;
   for (int i = 0; i < rows; i++)
   {
      cout << " " << output[i];
   }
   cout << endl;

   // ~~~~~~~~~~~~~~
   // Time for some polymorphism!
   // ~~~~~~~~~~~~~~
   // sparse_mat_poly is a pointer to a Matrix<double>
   // Because the CSRMatrix is a subclass of CSRMatrix that is ok
   Matrix<double> *sparse_mat_poly = new CSRMatrix<double>(rows, cols, nnzs, true);

   // We can even use sparse_mat_poly as inputs/outputs to routines that expect
   // What version of our matMatMult would be called here??
   // ~~~~~~~~
   // This line will segfault if called!
   // ~~~~~~~~   
   //sparse_mat_poly->matMatMult(*sparse_mat_poly, *sparse_mat_poly);

   // What if we did that on sparse_mat
   sparse_mat->matMatMult(*sparse_mat, *sparse_mat);

   // Can we do this??
   printStuffMatrix(*sparse_mat_poly);
   // Can we also do this??
   printStuffMatrix(*sparse_mat);

   delete sparse_mat_poly;
   delete sparse_mat;

	system("pause");

}