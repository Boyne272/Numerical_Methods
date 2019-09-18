/*
#include "header.h"
#include "Matrix.h"

void mat_prod(vector<double> &A, vector<double> &B, int row_A, int row_B, vector<double> &C);
void mat_trans(vector<double> &A, vector<double> &A_trans, int row_A);
void mat_inv(vector<double> &A, vector<double> &A_inv);
void linear_regress(vector<double> &data_x, vector<double> &data_y, double *a, double *b);


void test_mat_prod()
//int main()
// to test matrix multiplication
{
	// setup the non-square matrix to multiply
	vector<double> mat1 = { 0,1,2,3,4,5,6,7,8,9,10,11 };  // 4x3
	vector<double> mat2 = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 };  //3x5
	vector<double> res;

	// do the multiplication
	mat_prod(mat1, mat2, 4, 3, res);

	// print the result
	for (int i = 0; i < res.size(); i++)
	{
		cout << "element i = " << res[i] << endl;
	}
	system("pause");
}


void test_mat_inv()
// also test the product
{
	// setup the 2x2 matrix to invert
	vector<double> mat1 = { 1,2,3,4 };
	vector<double> res;
	vector<double> prod;

	// do inverse and product
	mat_inv(mat1, res);
	mat_prod(mat1, res, 2, 2, prod);

	// test product is identity as expected
	for (int i = 0; i < 4; i++)
	{
		if (i == 0 || i == 3)
		{
			assert(prod[i] == double(1));
		}
		else
		{
			assert(prod[i] == double(0));
		}
	}
	cout << "\ntest_mat_inv passed\n";
	system("pause");
}


void test_linear_reg()
{
	// initialise everything
	vector<double> vec_y = { 16,18,21,17,15, 12 };
	vector<double> vec_x = { 1,2,3,4,5,6 };
	double x_val;
	double result;
	double *res = &result;

	double a, b;

	//cout << "Give me a number bruh: ";
	//while (cin >> x_val) {
	//	lagrange(vec_x, vec_y, x_val, res);
	//	cout << "The value here is: " << result << endl << "Give me a numba bruh : ";
	//}

	linear_regress(vec_x, vec_y, &a, &b);

	system("pause");
}


int test_matrix_multiplication()
//int main()
{
	matrix A(3, 12), B(5, 15);

	for (int i = 0; i < 12; i++)
		A.append(i);
	for (int i = 0; i < 15; i++)
		B.append(i);

	A.print();
	B.print();

	matrix C = A * B;
	C.print();

	return 0;
}


int test_mat_inv()
//int main()
{
	// setup the 2x2 matrix to invert
	matrix mat(2, 4);
	for (int i = 0; i < 4; i++)
		mat.append(i*2);

	mat.print();
	mat.T().print();

	// do inverse and product
	matrix mat_inv = mat.inverse_2();
	matrix ident = mat_inv * mat;

	ident.print();

	system("pause");
	return 0;
}

*/