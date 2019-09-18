#include "header.h";
#include "matrix.h";
#include "menu.h";

	// I am aware that this is bad practice but steven did it in
	// lectures and without it I get a linker error
#include "generate_matrix.cpp";
#include "matrix.cpp";


	// menu for user interface
class matrix_menu : public menu {
private:
	string filename = "";
	ofstream save_file;

public:

		// these are all the commands for the unser interface
	matrix_menu() {
		this->command_map["example 0"] = 0;
		this->desctriptions["example 0"] = "create a simple dense matrix";
		this->command_map["example 1"] = 1;
		this->desctriptions["example 1"] = "create a simple sparse matrix (full)";
		this->command_map["example 2"] = 2;
		this->desctriptions["example 2"] = "create a simple sparse matrix (half-full)";
		this->command_map["example 3"] = 3;
		this->desctriptions["example 3"] = "example of dense matrix multiplication";
		this->command_map["example 3.5"] = 305;
		this->desctriptions["example 3.5"] = "example of sparse-dense matrix multiplication";
		this->command_map["example 4"] = 4;
		this->desctriptions["example 4"] = "example gauss-sidel solver for Ax=b with a 5x5 dense matrix";
		this->command_map["example 5"] = 5;
		this->desctriptions["example 5"] = "example gauss-sidel solver for Ax=b with a 12x12 sparse matrix";
		this->command_map["example 6"] = 6;
		this->desctriptions["example 6"] = "create positive definite sparse matrices of different sizes and same sparcity";
		this->command_map["example 7"] = 7;
		this->desctriptions["example 7"] = "create positive definite dense matrices of different sizes and same sparcity";
		this->command_map["time single"] = 8;
		this->desctriptions["time single"] = "timed gauss-sidel signle run on a random matrix for user specified conditions";
		this->command_map["time batch"] = 9;
		this->desctriptions["time batch"] = "timed gauss-sidel signle run on several random matrices with pre-determined conditions";
		this->command_map["open file"] = 10;
		this->desctriptions["open file"] = "create (or open pre-existing) a csv for saving batch solve results";
		this->command_map["close file"] = 11;
		this->desctriptions["close file"] = "close an open file";

		this->intro = "Matrix System Solver - ACSE 5.3 Coursework\nBy Richard Boyne CID 01057503";
	}


		// this sets up the problem then solves with gauss-seidel for a dense matrix
	void time_solve_dense(int iterations, int size, double sparcity,
					   	  int seed, int off_diag, bool print_mat,
						  bool save = 0) {
		
			// setup the problem
		Matrix_dense<double> A = rand_PDM_dense(size, seed, sparcity, off_diag);

			// the answer we want to find
		srand(seed);
		Matrix_dense<double> x(size, 1);
		for (int i = 0; i < size; i++)
			x.set(100 * double(rand()) / RAND_MAX);

			// the initial guess
		Matrix_dense<double> sol(size, 1, double(1));
			// the resultant vector for input
		Matrix_dense<double> B = A * x;

			// print the matrices
		if (print_mat) {
			cout << "\nMatrix A (dense)";
			A.print_matrix();
			A.info();
			cout << "\nMatrix B";
			B.print_matrix();
			cout << "\nMatrix guess";
			sol.print_matrix();
			cout << "\nMatrix x";
			x.print_matrix();
		}

			// solve
		cout << "solving dense ...\n";
		clock_t timer = clock();
		A.gauss_seidel(B, sol, iterations);
		timer = clock() - timer;
		cout << "ticks taken: " << timer << '\n';

			// find error
		double std = 0;
		for (int n = 0; n < size; n++)
			std += pow(x.i(n, 0) - sol.i(n, 0), 2);
		std = sqrt(std / size);
		cout << "standard deviation is: " << std << '\n';

			// save if wanted
		if (save) {
			this->save_file << "dense" << "," 
							<< size << ","
							<< sparcity << ","
							<< iterations << ","
							<< seed << ","
							<< std << ","
							<< timer << "\n";
			cout << "Saved to " << this->filename << '\n';
		}
	}


		// exactly the same as above but with a single line different to be sparse
		// can not for the life of me find a more elegent solution to just copying
		// the function and it eats me up inside knowing I have done this
	void time_solve_sparse(int iterations, int size, double sparcity,
				  		  int seed, int off_diag, bool print_mat,
						  bool save = 0) {

			// setup the problem (the different line)
		Matrix_sparse<double> A = rand_PDM_sparse(size, seed, sparcity, off_diag);

			// the answer we want to find
		srand(seed);
		Matrix_dense<double> x(size, 1);
		for (int i = 0; i < size; i++)
			x.set(100 * double(rand()) / RAND_MAX);

			// the initial guess
		Matrix_dense<double> sol(size, 1, double(1));
			// the resultant vector for input
		Matrix_dense<double> B = A * x;

			// print the matrices
		if (print_mat) {
			cout << "\nMatrix A (dense)";
			A.print_matrix();
			A.info();
			cout << "\nMatrix B";
			B.print_matrix();
			cout << "\nMatrix guess";
			sol.print_matrix();
			cout << "\nMatrix x";
			x.print_matrix();
		}

			// solve
		cout << "solving sparse ...\n";
		clock_t timer = clock();
		A.gauss_seidel(B, sol, iterations);
		timer = clock() - timer;
		cout << "ticks taken: " << timer << '\n';

			// find error
		double std = 0;
		for (int n = 0; n < size; n++)
			std += pow(x.i(n, 0) - sol.i(n, 0), 2);
		std = sqrt(std / size);
		cout << "standard deviation is: " << std << '\n';

			// save if wanted
		if (save) {
			this->save_file << "sparse" << ","
				<< size << ","
				<< sparcity << ","
				<< iterations << ","
				<< seed << ","
				<< std << ","
				<< timer << "\n";
			cout << "Saved to " << this->filename << '\n';
		}
	}


		// this is the code to be exicuted on each command
		// they are examples of how a user could use the code
	void commands(int opt) {
		switch(opt) {

			// example 0 - numberline dense matrix
		case 0: {
			Matrix_dense<double> example1(10, 10);
			for (int row = 0; row < example1.rows; row++)
				for (int col = 0; col < example1.cols; col++)
					example1.set(row, col, row * example1.cols + col);
			example1.print_matrix();
			example1.info();
			break;
		}


			// example 1 - numberline spares matrix
		case 1: {
			Matrix_sparse<double> example2(10, 10, 100);
			for (int row = 0; row < example2.rows; row++)
				for (int col = 0; col < example2.cols; col++)
					example2.set(row, col, row * example2.cols + col);
			example2.info();
			example2.print_matrix();
			break;
		}


			// example 2 - half filled sparse matrix
		case 2: {
			Matrix_sparse<double> example(10, 10, 50);
			for (int row = 0; row < example.rows; row++)
				for (int col = 0; col < example.cols; col++)
					if ((row + col) % 2 == 0)
						example.set(row, col, row * example.cols + col + 1);
			example.info();
			example.print_matrix();
			break;
		}


			// example 3 - dense*dense to dense matrix multiplication
		case 3: {
				// setup matrix A
			Matrix_dense<double> A(5, 5);
			for (int row = 0; row < A.rows; row++)
				for (int col = 0; col < A.cols; col++)
					A.set(row, col, row * A.cols + col);
			cout << "\nMatrix A:\n";
			A.print_matrix();

				// setup Identity matrix
			Matrix_dense<double> I(5, 5);
			for (int row = 0; row < I.rows; row++)
				for (int col = 0; col < I.cols; col++) {
					if (row == col)
						I.set(row, col, 1);
					else
						I.set(row, col, 0);
				}
			cout << "\nMatrix I:\n";
			I.print_matrix();

				// multiply A with identity matrix
			Matrix_dense<double> C = A * I;
			cout << "\nMatrix A*I\n";
			C.print_matrix();

				// multiply A by itself
			Matrix_dense<double> D = A * A;
			cout << "\nMatrix A*A\n";
			D.print_matrix();
			
			break;
		}

		
			// example 3.5 - sparse*dense to dense matrix multiplication
		case 305: {
				// setup matrix A
			Matrix_sparse<double> A(10, 10, 50);
			for (int row = 0; row < A.rows; row++)
				for (int col = 0; col < A.cols; col++)
					if ((row + col) % 2 == 0)
						A.set(row, col, row * A.cols + col + 1);
			cout << "\nMatrix A:\n";
			A.print_matrix();

				// setup Identity matrix
			Matrix_dense<double> I(10, 10);
			for (int row = 0; row < I.rows; row++)
				for (int col = 0; col < I.cols; col++) {
					if (row == col)
						I.set(row, col, 1);
					else
						I.set(row, col, 0);
				}
			cout << "\nMatrix I:\n";
			I.print_matrix();

				// multiply A with identity matrix
			Matrix_dense<double> C = A * I;
			cout << "\nMatrix A*I\n";
			C.print_matrix();

				// setup simple vector
			Matrix_dense<double> vec(10, 1);
			for (int row = 0; row < I.rows; row++)
				vec.set(row, 0, row);
			cout << "\nVector V:\n";
			vec.print_matrix();

				// multiply A with vec
			Matrix_dense<double> D = A * vec;
			cout << "\nMatrix A*A\n";
			D.print_matrix();

			break;
		}
				

			// example 4 - dense matrix gauss-seidel example
		case 4: {
				// setup positive definite matrix A
			Matrix_dense<double> A(5, 5);
			for (int row = 0; row < A.rows; row++)
				for (int col = 0; col < A.cols; col++) {
					if (col == row)
						A.set(row, col, row * A.cols + col + 20);
					else
						A.set(row, col, row * A.cols + col);
				}
			cout << "\nMatrix A:\n";
			A.print_matrix();

				// setup matrix B for which we know x = 2 everywhere
				// (found by doing matrix multiplication of Ax)
			Matrix_dense<double> B(5, 1);
			B.set(0, 0, 60.);
			B.set(1, 0, 110.);
			B.set(2, 0, 160.);
			B.set(3, 0, 210.);
			B.set(4, 0, 260.);
			cout << "\nMatrix B:\n";
			B.print_matrix();

				// setup an intial guess of x = 1 everywhere
			Matrix_dense<double> x(5, 1);
			for (int row = 0; row < x.rows; row++)
				x.set(row, 0, 1);
			cout << "\nIntial Guess:\n";
			x.print_matrix();

				// solve with gauss_seidel
			A.gauss_seidel(B, x, 100);
			cout << "\nResult after 100 iterations:\n";
			x.info();
			x.print_matrix();

			break;
		}


			// example 5 - sparse matrix gauss-seidel example
		case 5: {
				// setup positive definite matrix A with a reasonable sparsity
			Matrix_sparse<double> A(12, 12, 12 * 6);
			for (int row = 0; row < A.rows; row++)
				for (int col = 0; col < A.cols; col++) {
					if (col == row)  // make diagonals big
						A.set(row, col, row * A.cols + col + 15);
					else
						if ((row + col) % 2 == 0)  // make every other non-diagonal have value
							A.set(row, col, row * A.cols + col);
				}
			cout << "\nMatrix A:\n";
			A.print_matrix();

				// setup matrix B for which we know x = 2 everywhere
				// (found by doing matrix multiplication of Ax)
			Matrix_dense<double> B(12, 1);
			for (auto i : { 90., 246., 378., 534., 666., 822., 954., 1110., 1242., 1398., 1530., 1686. })
				B.set(i);
			cout << "\nMatrix B:\n";
			B.print_matrix();

				// setup an intial guess of x = 1 everywhere
			Matrix_dense<double> x(12, 1);
			for (int row = 0; row < x.rows; row++)
				x.set(row, 0, 1);
			cout << "\nIntial Guess:\n";
			x.print_matrix();

				// solve with gauss_seidel
			A.gauss_seidel(B, x, 100);
			cout << "\nResult after 100 iterations:\n";
			x.info();
			x.print_matrix();

			break;
		}
		

			// example 6 - Generation of random sparse matices with fixed sparcity
		case 6: {
			for (int size : {10, 20, 30, 40}) {
				Matrix_sparse<double> A = rand_PDM_sparse(size, 20, 0.6);
				A.print_matrix();
				A.info();
			}
			break;
		}

	
			// example 7 - Generation of random dense matices with fixed sparcity
		case 7: {
			for (int size : {10, 20, 30, 40}) {
				Matrix_dense<double> A = rand_PDM_dense(size, 20, 0.6);
				A.print_matrix();
				A.info();
			}
			break;
		}
		

			// generate and time either dense or sparse matrices with user given parameters
		case 8: {

				// vairables (to be user inputs)
			cout << "Number of iterations: ";
			int iterations = this->enter_int(1, 1000);

			cout << "Matrix size: ";
			int size = this->enter_int(1, 10000);

			cout << "Matrix sparcity (1-100): ";
			const int min_sparse = int(100 / size);
			double sparcity = double(this->enter_int(min_sparse, 100))/100;

			cout << "Random generator seed: ";
			int seed = this->enter_int();

			cout << "Max off diagonal value: ";
			int off_diag = this->enter_int(0, 1000);

			cout << "Print matrices (y/n): ";
			bool print_mat = this->enter_bool();

			cout << "Save results (ensure file is opened)? (y/n): ";
			bool save = this->enter_bool();

			cout << "Sparse, Dense or both matrix type (s/d/b): ";
			char opts[] = { 's', 'd', 'b'};
			char mat_type = this->enter_options<char>(opts, 3);

			if (mat_type == 'd' || mat_type == 'b')
				this->time_solve_dense(iterations, size, sparcity, seed,
									   off_diag, print_mat, save);
			if (mat_type == 's' || mat_type == 'b')
				this->time_solve_sparse(iterations, size, sparcity, seed,
									    off_diag, print_mat, save);

			break;
		}


			// generate and time a large collection of matrices
		case 9: {
				// check a save file exists
			if (!this->save_file.is_open()) {
				cout << "Error: A file to save in must be opened with open file\n";
				break;
			}

				// set ip the parameters to generate and run with
			const int iterations = 100;
			const int off_diag = 5;
			int sizes[10];
			for (int i = 0; i < 10; i++)
				sizes[i] = 10 * pow(2, i);

				// get user choices
			cout << "Matrix sparcity (1-100): ";
			const int min_sparse = int(100 / sizes[0]);
			double sparcity = double(this->enter_int(min_sparse, 100)) / 100;

			cout << "Number of repeats: ";
			const int num_repeats = this->enter_int(0, 100);
			int* seeds = new int[num_repeats];
			for (int i = 0; i < num_repeats; i++)
				seeds[i] = i;

				// print what is beeing done
			cout << "Will solved for " << num_repeats << " repeats with sparcity " << sparcity << " for sizes: ";
			for (int i = 0; i < 10; i++)
				cout << sizes[i] << " ";
			cout << '\n';
			system("pause");

				// solve for each and save results
			for (auto size : sizes)
				for (int i=0; i < num_repeats; i++) {
					cout << "\nSize " << size << " seed " << seeds[i] << '\n';
					this->time_solve_dense(iterations, size, sparcity, seeds[i],
						off_diag, 0, 1);
					this->time_solve_sparse(iterations, size, sparcity, seeds[i],
						off_diag, 0, 1);
				}
			cout << "finished solving run\n";

			break;
		}


			// open save file
		case 10: {
			cout << "Give a file name: ";
			this->filename = this->enter_string(".csv", " ");

			this->save_file.open(this->filename, std::ios_base::app);

			if (this->save_file.is_open())
				cout << this->filename << " successfully opened\n";
			else
				cout << "Error: could not open file\n";

				// set the heading for saving data
			cout << "Would you like to create a header (y/n): ";
			if (this->enter_bool()) {
				this->save_file << "Type,Size,Density,Iterations,Seed,STD,Time\n";
				cout << "File header created\n";
			}
			break;
		}


			// close save file
		case 11: {
				// if the file is opened try to close it
			if (this->save_file.is_open()) {
				this->save_file.close();

					// if it was not closed give an error message
				if (!this->save_file.is_open()) {
					cout << this->filename << " successfully closed\n";
					this->filename = "";
				}
				else
					cout << "Error: could not close file";
			}
			else
				cout << "Error: No file is currently open";
			break;
		}

		}
	}
};


int main() {
		// create and start the uwer interface
	matrix_menu inst = matrix_menu();
	inst.start();
}