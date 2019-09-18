#include "header.h"
#include "Matrix.h"

double error_overall(matrix &observed, matrix &forecast)
// matrix &observed = temperature and years of input file eg JMA_tokyo_data.csv (at specified month)
// matrix &forecast = temperature and years within model range (at specified month)
{
	//assert(observed.size == forecast.size);
	double standard_error = 0;

	for (int i = 0; i < observed.n_row(); i++)
	{
		//cout << "M[" << i << "] = " << observed.i(i, 0) << "\n";
		standard_error += pow((observed.i(i, 1) - forecast.i(i, 1)), 2);
	}
	standard_error = standard_error/ observed.n_row();
	return sqrt(standard_error);
}


void error_per_point(matrix &pred, matrix &orig,
					 matrix &error, double model_data_end)
// matrix &pred = temperature and years within model range (at specified month)
// matrix &orig = temperature and years of input file eg JMA_tokyo_data.csv (at specified month)
// matrix &error = matrix to store errors
// double model_data_end = prediction year as inputted by user
{
	// go back through original data set, and compare with
	// data outside of specified model
	// eg model is for 1900 - 2010, prediction wanted is 2050
	// compare predicted values outside of 2010 to original data
	assert(error.size == 0);
	error.n_col = 2;

	for (int i = 0; i < orig.n_row(); i++)
	{
		// find the corresponding fitting value
		for (int j = 0; j < pred.n_row(); j++)
		{
			if (pred.i(j, 0) == orig.i(i, 0) && pred.i(j, 0) > model_data_end)
			{
				error.append(pred.i(j, 0));
				error.append(abs(pred.i(j, 1) - orig.i(i, 1)));
				break;
			}
		}
	}
}

