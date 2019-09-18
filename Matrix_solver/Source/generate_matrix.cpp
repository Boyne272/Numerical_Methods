#include "header.h";
#include "matrix.h";
#include "menu.h";


// function to generate a random positive (dense) matrix with a given sparcity
Matrix_dense<double> rand_PDM_dense(int size, int seed = 10, double sparcity = 1.00,
	int max_off_diag = 3)
	// Here I use knowledge that positive diagonally dominant matrices are
	// guarenteed positive definite and fill values at even intervals to
	// achive as close to desired sparcity as possible
{

	// validate inputs 
	assert(sparcity >= (1 / size));
	assert(sparcity <= 1);

	// set seed
	srand(seed);
	// useful constants
	const int n = size * size * sparcity;  // number of elements to set
	const int non_diag = n - size;  // number of non-diagonal elements

		// setup the matrix and initialise to zero everywhere
	Matrix_dense<double> A(size, size, double(0));

	// find and set the off diagonal values to use
	double step = double((size*size) - size) / non_diag;  // the average distance between two off_diagonals
	int * off_diags = new int[non_diag];
	int index;
	for (int i = 0; i < non_diag; i++) {
		off_diags[i] = (rand() % max_off_diag) + 1;
		index = step * i;
		index += int(index / size) + 1;  // this finds the off-diagonal entry to use
		const int row = int(index / size);
		const int col = index % size;
		A.set(row, col, off_diags[i]);
	}

	// find the diagonal values
	const double row_length = non_diag / size;
	double val;
	for (int row = 0; row < size; row++) {
			// give it a non-zero value incase row is empty
		val = (rand() % max_off_diag) + 1;
		for (int j = row * row_length; j < (row + 1) * row_length; j++)
			val += off_diags[j];
		A.set(row, row, val);
	}

	// delete working arrays
	delete[] off_diags;

	return A;
}


// function to generate a random positive (sparse) matrix with a given sparcity
Matrix_sparse<double> rand_PDM_sparse(int size, int seed = 10, double sparcity = 1.00,
	int max_off_diag = 3)
	// Here I use knowledge that positive diagonally dominant matrices are
	// guarenteed positive definite and fill values at even intervals to
	// achive as close to desired sparcity as possible
{

	// validate inputs 
	assert(sparcity > (1 / size));
	assert(sparcity <= 1);

	// set seed
	srand(seed);
	// useful constants
	const int n = size * size * sparcity;  // number of elements to set
	const int non_diag = n - size;  // number of non-diagonal elements

		// setup the matrix
	Matrix_sparse<double> A(size, size, n);

	// find the off diagonal values to use
	int * off_diags = new int[non_diag];
	for (int i = 0; i < non_diag; i++) {
		off_diags[i] = (rand() % max_off_diag) + 1;
	}

	// find the diagonal values
	const double row_length = non_diag / size;
	int * diag = new int[size];
	for (int row = 0; row < size; row++) {
			// give it a non-zero value incase row is empty
		diag[row] = (rand() % max_off_diag) + 1;
		for (int j = row * row_length; j < (row + 1) * row_length; j++)
			diag[row] += off_diags[j];
	}

	// set the matrix values
	// this gets quite complex as CSR format requiers we enter values
	// in order of occurance in the matrix, hence we need to input the
	// diagonals into the matrix.set method in at the right time (not easy)
	int index;
	int act_index;
	int next_diag = 0;
	int count_diag = 0;
	double step = double((size*size) - size) / non_diag;  // the average distance between two off_diagonals

	for (int i = 0; i < non_diag; i++) {
		index = step * i;
		// I found next line by playing in excel till it worked
		// it add the correct amount to curr_index so that it
		// skips all diags up to this index
		act_index = index + int(index / size) + 1;

		// if we crossed a diagonal entery the diagonal first 
		while (act_index > next_diag) {
			A.set(count_diag, count_diag, diag[count_diag]);
			count_diag++;
			next_diag += size + 1;
		}

		const int row = int(act_index / size);
		const int col = act_index % size;
		A.set(row, col, off_diags[i]);
	}

	// set any diagonals we never crossed
	while (next_diag < (size*size)) {
		A.set(count_diag, count_diag, diag[count_diag]);
		count_diag++;
		next_diag += size + 1;
	}

	// delete working arrays
	delete[] off_diags;
	delete[] diag;

	return A;
}
