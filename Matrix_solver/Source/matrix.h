#pragma once

// Top quote: "if a mistake is made you want it to break and break hard" - Steven

template <class T>
class Matrix
{
public:
		// initialise vairables
	int rows = -1;
	int cols = -1;
	int size = -1;


		// default constructor
	Matrix();

		// assignment of defualt values for intialisation list
	Matrix(int rows, int cols);
	
		// default destructor
	~Matrix();

		// requiered functions
	virtual T i(int row, int col) = 0;  // index an element
	virtual void set(int row, int col, T val) = 0;  // set ith element

		// print readable format
	virtual void print_matrix();

		// print matrix info
	virtual void info();
};


template <class T>
class Matrix_dense : public Matrix<T>
{
public:
		// data storage method
	unique_ptr<T[]> uptr_data; // handles memory
	T* ptr_data = nullptr;  // indexable version
	int entry_counter = 0;


		// contructor: memory assigned by the constructor
	Matrix_dense(int rows, int cols);

		// contructor: memory assigned by the constructor with initial value set
	Matrix_dense(int rows, int cols, T init_val);

		// constructor: memory pre-assigned (we wont handle the memory)
	Matrix_dense(int rows, int cols, T* value_ptr);


		// index the ith element
	T i(int row, int col);


		// set ith element
	void set(int row, int col, T val);


	// set the next element for sequential initialisation 
	void set(T val);


		// dense matrix multiplication with dense matrix to dense matrix
	Matrix_dense<T> operator*(Matrix_dense<T> &other_mat);


		// solve Ax=b with gauss seidel method
	void gauss_seidel(Matrix_dense &b, Matrix_dense &guess, int iterations);		

		// print matrix info
	void info();

};


template <class T>
class Matrix_sparse : public Matrix<T>
{
public:
		// unique pointers to handle memory
	unique_ptr<T[]> uptr_value;
	unique_ptr<int[]> uptr_col_index;
	unique_ptr<int[]> uptr_row_pos;

		// raw pointers for indexing
	T* ptr_value = nullptr;
	int* ptr_col_index = nullptr;
	int* ptr_row_pos = nullptr;

		// non-zero matrix size
	int nz_size = -1;

		// setting up matrix counters
	int entery_counter = 0;
	int row_counter = 0;


		// contructor: memory assigned by the constructor
	Matrix_sparse(int rows, int cols, int number_non_zero);

		// index the ith element
	T i(int row, int col);
	// this method allows for backwards compatibility with any dense matrix method
	// however it is very inefficent and should not be used for optimised algs

		// set ith element
	void set(int row, int col, T val);
	// Data must be set in order of occurance in a row stored matrix


		// solve Ax=b with gauss seidel method - optimised for sparse matrix
	void gauss_seidel(Matrix_dense<T> &b, Matrix_dense<T> &guess, int iterations);
	// ** matrix assumed positive definite with non-zero diagonals

		// sparse matrix multiplication with dense matrix to dense matrix
	Matrix_dense<T> operator*(Matrix_dense<T> &other_mat);


		// print readable format with * for zero enteries
	void print_matrix(); 


		// print matrix info
	void info();

};


// below is supposed to be the declaration of rand_PDM* functions but
// every time I try to do this I get a linker error :'(

	// function to generate a random positive (dense) matrix with a given sparcity
	//Matrix_dense<double> rand_PDM_dense(int size, int seed = 10, double sparcity = 1.00,
	//	int max_off_diag = 3);
		// Here I use knowledge that positive diagonally dominant matrices are
		// guarenteed positive definite and fill values at even intervals to
		// achive as close to desired sparcity as possible


	// function to generate a random positive (sparse) matrix with a given sparcity
	//Matrix_sparse<double> rand_PDM_sparse(int size, int seed = 10, double sparcity = 1.00,
	//	int max_off_diag = 3);
		// Here I use knowledge that positive diagonally dominant matrices are
		// guarenteed positive definite and fill values at even intervals to
		// achive as close to desired sparcity as possible