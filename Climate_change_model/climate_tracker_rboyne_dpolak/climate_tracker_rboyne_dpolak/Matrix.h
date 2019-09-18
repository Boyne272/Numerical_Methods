#pragma once

using namespace std;

class matrix
	// only allows a 2d matrix
{
public:

	// declare attributes
	vector<double> vec;
	int n_col, size;
	bool filled = false;

	matrix(int cols = 1, int init_size = 10)
	{
		vec.reserve(init_size);
		n_col = cols;
		size = 0;
	}

	~matrix()
	{
		//cout << "\nmatrix deconstructed\n";
		vec.clear();
	}

	void append(double val)
	// add a value and check if all rows are full
	// double value = value which to append
	{
		vec.push_back(val);
		size++;
		if (vec.size() % n_col == 0)
			filled = true;
		else
			filled = false;
	}

	int n_row()
		// find how many rows currently in the matrix
	{
		return int(vec.size() / n_col);
	}

	double i(int i, int j)
	// get the n, m element
	// int i = number of rows
	//int j = number of columns
	{
		return vec[i*n_col + j];
	}

	void set(int i, int j, double val)
	// int i = number of rows
	//int j = number of columns
	// double value = value which to equal 
	{
		vec[i*n_col + j] = val;
	}

	double from_end(int n = 0, int col = 0)
	// index from the end of the matrix
	// int n = row number
	// col = 0 -> year, 1 -> temperature
	{ 
		return i(n_row()-1-n, col);
	}

	void print()
	// pretty display of the matrix
	{
		cout << "Matrix (" << n_row() << ", " << n_col << ")\n";
		for (int i = 0; i < n_row(); i++)
		{
			for (int j = 0; j < n_col; j++)
				cout << "\t" << vec[i*n_col + j] << "\t";
			cout << endl;
		}
	}

	matrix operator*(matrix& other_matrix)
		// matrix multiplication
		// matrix& other_matrix = matrix to be multiplied
	{
		// can't multiply unfilled matrix's
		assert(filled);
		assert(other_matrix.filled);

		// must have correct dimensions
		assert(n_col == other_matrix.n_row());

		// create result matrix
		matrix result(other_matrix.n_col, n_row()*other_matrix.n_col);

		// declare tmp variables for the loop
		double tmp;

		//loop over every element in the output vector
		for (int i = 0; i < n_row(); i++)
		{
			for (int j = 0; j < other_matrix.n_col; j++)
			{
				tmp = 0;
				// loop over every element in the sum
				for (int n = 0; n < n_col; n++)
				{
					tmp += vec[i*n_col + n] * other_matrix.vec[n*other_matrix.n_col + j];
				}
				// add this to the final array
				result.append(tmp);
			}
		}

		// check the final matrix makes sense
		assert(result.filled);

		return result;
	}

	matrix T()
		// find the transpose
	{
		assert(filled);
		// define local vairables
		double store;
		matrix trans(n_row(), vec.size());

		// loop over every element
		for (int i = 0; i < n_col; i++)
		{
			for (int j = 0; j < n_row(); j++)
			{
				// add the corresponding element
				trans.append(vec[j*n_col + i]);
			}
		}

		return trans;
	}

	matrix inverse_2()
		// inverse a 2x2 matrix
	{
		assert(vec.size() == 4);
		assert(n_col == 2);

		// setup
		double det = (vec[0] * vec[3] - vec[1] * vec[2]);
		matrix inv(2, 4);

		// for each element enter its value
		inv.append(vec[3] / det);
		inv.append(-vec[1] / det);
		inv.append(-vec[2] / det);
		inv.append(vec[0] / det);

		return inv;
	}
};