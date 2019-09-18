#include "header.h"
#include "Matrix.h"


double lagrange(matrix &data, double x)
// matrix &data = temperature and years within model range (at specified month)
// double x = point to evaluate at (prediction value)
{
	// set up
	double tmp, result(0);

	// validate input
	assert(data.filled);

	// for every data set pair
	for (int i = 0; i < data.n_row(); i++)
	{
		tmp = 1;
		for (int j = 0; j < data.n_row(); j++)
		{
			if (i != j)
				// multiply the term to the product
				tmp = tmp * (x - data.i(j, 0)) / (data.i(i, 0) - data.i(j, 0));
		}

		// add to the result
		result = result + tmp * data.i(i, 1);
	}

	return result;
}


void linear_regress(matrix &data, double *a, double *b)
// find the linear_regression (i.e. least square fit) 
// by matrix method taught in lect 2, pass in the x and y data
// matrix &data = temperature and years within model range (at specified month)
// double a, double b = variables in which to store fitted variables
{
	//validate input
	assert(data.filled);

	// setup G and d matrixs
	matrix G(2, data.size);
	matrix d(1, data.n_row());
	for (int i = 0; i < data.n_row(); i++)
	{
		G.append(double(1));
		G.append(data.i(i, 0));
		d.append(data.i(i, 1));
	}

	// solve the equation
	matrix tmp1, tmp2, model;
	//vector<double> G_T, G_prod, G_prod_inv, G_T_data, model;

	// find inverse square mat
	tmp1 = G.T() * G;
	tmp1 = tmp1.inverse_2();

	// find G_T data transpose
	tmp2 = G.T() * d;

	// find the model parameters
	model = tmp1 * tmp2;

	*a = model.i(0, 0);
	*b = model.i(0, 1);
	cout << "\nFitted Values:\n" << "\tIntercept = " << *a
		 << "\n\tGradient = " << *b << endl;
}

