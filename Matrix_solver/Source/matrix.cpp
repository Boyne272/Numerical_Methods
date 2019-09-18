#include "header.h"
#include "matrix.h"

// -------------------------------------------------------
//		Matrix Class Functions
// -------------------------------------------------------

	// default constructor
template <class T>
Matrix<T>::Matrix() {
	//cout << "constructor at " << this << " called (default for Matrix)\n";
}


	// assignment of defualt values for intialisation list
template <class T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols), size(rows*cols) {
	//cout << "constructor at " << this << " called (non-default for Matrix)\n";
}


	// default destructor
template <class T>
Matrix<T>::~Matrix() {
	//cout << "destructor  at " << this << " called (default for Matrix) \n";
}


	// print readable format
template <class T>
void Matrix<T>::print_matrix() {
	cout << "\nMatrix Values\n-------------\n";
	cout << "Dimensions (" << this->rows << ", " << this->cols << ")\n";
	for (int row = 0; row < this->rows; row++) {
		for (int col = 0; col < this->cols; col++)
			cout << "\t" << this->i(row, col);
		cout << "\n\n";
	}
}


	// print matrix info
template <class T>
void Matrix<T>::info() {
	cout << "\nMatrix Info\n"
		<< "-----------\n"
		<< "Default Matrix\n"
		<< "Location " << this << '\n'
		<< "Dimensions (" << this->rows << ", " << this->cols << ")\n";
}



// -------------------------------------------------------
//		Matric_dense Class Functions
// -------------------------------------------------------


	// contructor: memory assigned by the constructor
template <class T>
Matrix_dense<T>::Matrix_dense(int rows, int cols) : Matrix<T>(rows, cols) {
	//cout << "constructor at " << this << " called (memory assigning for Matrix_dense)\n";
		// set up the pointers
	uptr_data.reset(new T[this->size]);
	ptr_data = uptr_data.get();
}


	// contructor: memory assigned by the constructor with initial value set
template <class T>
Matrix_dense<T>::Matrix_dense(int rows, int cols, T init_val) : Matrix<T>(rows, cols) {
	//cout << "constructor at " << this << " called (default value for Matrix_dense)\n";
		
		// set up the pointers
	uptr_data.reset(new T[this->size]);
	ptr_data = uptr_data.get();

		// set all intial values to zero
	for (int i = 0; i < this->size; i++)
		ptr_data[i] = init_val;
}


	// constructor: memory pre-assigned (we wont handle the memory)
template <class T>
Matrix_dense<T>::Matrix_dense(int rows, int cols, T* value_ptr) : Matrix<T>(rows, cols) {
	//cout << "constructor at " << this << " called (pre-assigned for Matrix_dense)\n";
	ptr_data = value_ptr;
	this->preallocated = true;
}


	// index the ith element
template <class T>
T Matrix_dense<T>::i(int row, int col) {
	return this->ptr_data[row * this->cols + col];
}


	// set ith element
template <class T>
void Matrix_dense<T>::set(int row, int col, T val) {
	this->ptr_data[row * this->cols + col] = val;
}


	// set the next element for sequential initialisation
template <class T>
void Matrix_dense<T>::set(T val) {
	this->ptr_data[this->entry_counter] = val;
	this->entry_counter++;
}


	// dense matrix multiplication with dense matrix to dense matrix
template <class T>
Matrix_dense<T> Matrix_dense<T>::operator*(Matrix_dense<T> &other_mat)
{
	// check matrix dimensions are correct
	if (this->cols != other_mat.rows)
		cerr << "Matrix dimension mismatch ("
		<< this->rows << ", " << this->cols << " with "
		<< other_mat.rows << ", " << other_mat.cols << "\n\n";

	// create matrix to store result
	Matrix_dense<T> result(this->rows, other_mat.cols);
	T tmp;

	// loop over every elements in result
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			tmp = 0;
			// loop over the sum for the element
			for (int n = 0; n < this->cols; n++)
				tmp += this->i(i, n) * other_mat.i(n, j);
			result.set(i, j, tmp);
		}
	}
	return result;
}


	// solve Ax=b with gauss seidel method
template <class T>
void Matrix_dense<T>::gauss_seidel(Matrix_dense &b, Matrix_dense &guess, int iterations)
	// ** matrix assumed positive definite with non-zero diagonals
{
	// check b and guess is of correct dimensions
	if (b.cols != 1)
		cerr << "\n***Matrix b is not a column vector***\n\n";
	if (guess.cols != 1)
		cerr << "\n***Matrix guess is not a column vector***\n\n";
	if (guess.rows != b.rows)
		cerr << "\n***Matrix guess has different number of rows to matrix b***\n\n";
	if (b.rows != this->cols)
		cerr << "\n***Matrix dimension mismatch ("
		<< this->rows << ", " << this->cols << " with "
		<< b.rows << ", " << b.cols << "***\n\n";

	// setup some vairables used in iterations
	const int n = b.size;
	T* res = guess.ptr_data;
	double tmp;

	// loop over iterations
	for (int iter = 0; iter < iterations; iter++) {
		// loop over every value in result
		for (int i = 0; i < n; i++) {
			const int row_index = i * this->cols;
			tmp = b.i(i, 0);
			// loop over every element in the row of this
			for (int j = 0; j < n; j++) {
				if (i != j) {
					// optimisation note - the below line is slower
					// as compiler can't optimise
					//tmp -= this->i(i, j) * res[j];
					tmp -= this->ptr_data[row_index + j] * res[j];
				}
			}
			//res[i] = tmp / this->i(i, i);
			res[i] = tmp / this->ptr_data[row_index + i];
		}
	}

	return;
}


	// print matrix info
template <class T>
void Matrix_dense<T>::info() {
	cout << "\nMatrix Info\n"
		<< "-----------\n"
		<< "Dense Matrix\n"
		<< "Location " << this << '\n'
		<< "Dimensions (" << this->rows << ", " << this->cols << ")\n";
	int nzero_cout = 0;
	for (int i = 0; i < this->size; i++)
		if (this->ptr_data[i] != T(0))
			nzero_cout++;
	cout << "Non-zero enteries " << nzero_cout << " ("
		<< 100 * nzero_cout / double(this->size) << "%)\n\n";
}



// -------------------------------------------------------
//		Matric_sparse Class Functions
// -------------------------------------------------------


	// contructor: memory assigned by the constructor
template <class T>
Matrix_sparse<T>::Matrix_sparse(int rows, int cols, int number_non_zero) :
	Matrix<T>(rows, cols), nz_size(number_non_zero) {
	//cout << "constructor at " << this << " called (memory assigning for Matrix_sparse)\n";

		// set the unique pointers
	uptr_value.reset(new T[this->nz_size]);
	uptr_col_index.reset(new int[this->nz_size]);
	uptr_row_pos.reset(new int[this->rows + 1]);

	// set the raw pointers
	ptr_value = uptr_value.get();
	ptr_col_index = uptr_col_index.get();
	ptr_row_pos = uptr_row_pos.get();

	// initialise the rows array for the set method to work
	ptr_row_pos[0] = 0;
	for (int i = 1; i < this->rows + 1; i++)
		ptr_row_pos[i] = this->nz_size;
}


// index the ith element
template <class T>
T Matrix_sparse<T>::i(int row, int col)
	// this method allows for backwards compatibility with any dense matrix method
	// however it is very inefficent and should not be used for optimised algs
{
	// find the element indexs of this row
	const int start = this->ptr_row_pos[row];
	const int end = this->ptr_row_pos[row + 1];

	// loop over every element in the row
	for (int i = start; i < end; i++)
		if (this->ptr_col_index[i] == col)
			return this->ptr_value[i];

	return 0;
}


// set ith element
template <class T>
void Matrix_sparse<T>::set(int row, int col, T val)
	// Data must be set in order of occurance in a row stored matrix
{
	this->ptr_value[entery_counter] = val;
	this->ptr_col_index[entery_counter] = col;

	while (this->row_counter < row)
	{
		this->row_counter++;
		this->ptr_row_pos[this->row_counter] = this->entery_counter;
	}

	this->entery_counter++;

	return;
	//this->data_ptr[row * this->cols + col] = val;
}


// solve Ax=b with gauss seidel method - optimised for sparse matrix
template <class T>
void Matrix_sparse<T>::gauss_seidel(Matrix_dense<T> &b, Matrix_dense<T> &guess, int iterations)
	// ** matrix assumed positive definite with non-zero diagonals
{
	// check b and guess is of correct dimensions
	if (b.cols != 1)
		cerr << "\n***Matrix b is not a column vector***\n\n";
	if (guess.cols != 1)
		cerr << "\n***Matrix guess is not a column vector***\n\n";
	if (guess.rows != b.rows)
		cerr << "\n***Matrix guess has different number of rows to matrix b***\n\n";
	if (b.rows != this->cols)
		cerr << "\n***Matrix dimension mismatch ("
		<< this->rows << ", " << this->cols << " with "
		<< b.rows << ", " << b.cols << "***\n\n";

	// setup some vairables used in iterations
	const int n = b.size;
	T* res = guess.ptr_data;
	double tmp;
	T diag_element;

	// loop over iterations
	for (int iter = 0; iter < iterations; iter++) {
		// loop over every value in result
		for (int i = 0; i < n; i++) {
			tmp = b.i(i, 0);
			// loop over every element in the row of this
			// here this is done by finding the indexs in values that correspond
			// to row i, then setting the column j to the column recoreded
			// for that entry and using it in the loop as before
			const int start = this->ptr_row_pos[i];
			const int end = this->ptr_row_pos[i + 1];
			for (int index = start; index < end; index++) {
				const int j = this->ptr_col_index[index];
				if (i != j)
					tmp -= this->ptr_value[index] * res[j];
				else
					diag_element = this->ptr_value[index];  // stored for the next line
			}
			res[i] = tmp / diag_element;
		}
	}

	return;
}


// sparse matrix multiplication with dense matrix to dense matrix
template <class T>
Matrix_dense<T> Matrix_sparse<T>::operator*(Matrix_dense<T> &other_mat)
{
	// check matrix dimensions are correct
	if (this->cols != other_mat.rows)
		cerr << "Matrix dimension mismatch ("
		<< this->rows << ", " << this->cols << " with "
		<< other_mat.rows << ", " << other_mat.cols << "\n\n";

	// create matrix to store result with zeros in intiially
	Matrix_dense<T> result(this->rows, other_mat.cols, T(0));
	T tmp;

	// loop over every elements in result
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			tmp = 0;
			// loop over the sum for the element
			// only need consider elements of sum where sparse is non-zero
			const int start = this->ptr_row_pos[i];
			const int end = this->ptr_row_pos[i + 1];
			for (int index = start; index < end; index++) {
				const int n = this->ptr_col_index[index];
				tmp += this->ptr_value[index] * other_mat.i(n, j);
			}
			result.set(i, j, tmp);
		}
	}
	return result;
}


// print readable format with * for zero enteries
template <class T>
void Matrix_sparse<T>::print_matrix() {
	cout << "\nMatrix Values\n-------------\n";
	cout << "Dimensions (" << this->rows << ", " << this->cols << ")\n";
	for (int row = 0; row < this->rows; row++) {
		for (int col = 0; col < this->cols; col++) {
			const T val = this->i(row, col);
			if (val != T(0))
				cout << "\t" << this->i(row, col);
			else
				cout << "\t" << "*";
		}
		cout << "\n\n";
	}
}


// print matrix info
template <class T>
void Matrix_sparse<T>::info() {
	cout << "\nMatrix Info\n"
		<< "-----------\n"
		<< "Spares Matrix\n"
		<< "Location " << this << '\n'
		<< "Dimensions (" << this->rows << ", " << this->cols << ")\n"
		<< "Non-zero enteries " << this->nz_size << " (" << this->nz_size / double(this->size) * 100 << "%)\n"
		<< "Value Array: \n";
	for (int i = 0; i < this->nz_size; i++)
		cout << " " << this->ptr_value[i];
	cout << "\nColumn Array: \n";
	for (int i = 0; i < this->nz_size; i++)
		cout << " " << this->ptr_col_index[i];
	cout << "\nRow Array: \n";
	for (int i = 0; i < this->rows + 1; i++)
		cout << " " << this->ptr_row_pos[i];
	cout << "\n\n";
}
