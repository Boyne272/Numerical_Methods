#include "header.h"
#include "Matrix.h"


//
// Declare functions
//
// load data
void load_data(string filepath, int column, matrix &mat);
void truncate_data(matrix &data_trunc, double &start_year, double &end_year, matrix &data_orig);
//line fitting
void linear_regress(matrix &data, double *a, double *b);
double lagrange(matrix &data, double x);
void init_month(map<const int, string> &map_month);
// Plotting
void plot(string plot_string, string title = "Tokyo", string xlabel = "Time (years)", string ylabel = "Temp (Deg.C)");
string save_to_csv(matrix &mat, string prefix = "data_", string header = "");
//Choice Four and Three Erro fitting
double error_overall(matrix &observed, matrix &forecast);
void error_per_point(matrix &pred, matrix &orig, matrix &error, double model_data_end);


// input validation
template<typename T>
T input(T max, T min)
/* Validates the user input depending on type and range*/
{
	T value;
	while (true)
	{
		// get user input
		cin >> value;

		// if the input method failed (eg. wrong data type)
		if (!(cin))
		{
			cout << "\nWrong input type, please try again: ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
		}
		// if the value is out of range
		else if (value < min || value > max)
		{
			cout << "\nValue is out of range, please try again: ";
		}
		// else return the input
		else
		{
			return value;
		}
	}
}


void choice_one(matrix &data, int &pred_month)
// matrix &data = temperature and years within model range (at specified month)
// int &pred_month = prediction month as inputed by user
{
	// find linear regression parameters
	double a, b;
	linear_regress(data, &a, &b);

	// setup vairables
	double result_lin, result_lag, pred_year;
	char quit;
	cout << "\n\n Wanted prediction year (capped at 2100) : ";
	while (true)
	{
		pred_year = input<double>(2100, data.i(0, 0));
		// find the values at pred_year
		result_lag = lagrange(data, pred_year);
		result_lin = a + b * pred_year;

		// output answer
		cout << "\nPredictions for " << pred_month << "/" << pred_year << endl;
		cout << "\tLinear Regression: " << result_lin << " degree celsius average" << endl;
		cout << "\tLagrange Interpolation: " << result_lag << " degree celsius average" << endl;
		
		// repeat
		cout << "Enter q to quit or any other character predict again: ";
		cin >> quit;
		if (quit == 'q')
			return;
		else
			cout << "\n\nPrediction year: ";
	}
}

string choice_two(matrix &data, double pred_year, string name_prefix = "data_")
// matrix &data = temperature and years within model range (at specified month)
// int &pred_year = prediction year 
//  string name_prefix = (optional) .csv file name prefix
{
	//find linear regression parameters
	double a, b;
	linear_regress(data, &a, &b);

	// create vector of predictions
	matrix pred(2);
	for (double i = data.i(0, 0); i <= pred_year; i++)
	{
		// add the year
		pred.append(i);
		// add the prediction temp
		pred.append(a + b * i);
	}

	// save and return the file_name
	return save_to_csv(pred, name_prefix);
}


void choice_three(matrix &data, int month, string original_file,
				  double pred_year, bool keep_save = 0)
// matrix &data = temperature and years within model range (at specified month)
// int month = prediction month as inputed by user
// string original file = eg JMA_tokyo_data.csv, to plot raw data with prediction data
// bool keep_save = (optional) if true, .csv file created are saved
{
	// create the fitting and get the file name
	string file_name = choice_two(data, pred_year);

	// plot it
	string plot_command = "plot '" + file_name + "' using 1:2 with lines title 'Regressions', '" +
						  //"[" + to_string(int(data.i(0, 0))) + "," + to_string(int(pred_year)) + "] '" +
						  original_file + "' every ::1 using 1:" +
						  to_string(month + 1) + " title 'Data'\n";
	plot(plot_command, "Tokyo Linear Regression");

	// remove file if not wanted
	if (!keep_save)
	{
		cout << "removing file " << file_name << endl;
		remove(file_name.c_str());
	}
 }

// workout the overall error
void choice_four(matrix &data, string original_file, double final_pred_year, int month)
// matrix &data = temperature and years within model range (at specified month)
// string original file = eg JMA_tokyo_data.csv, to plot raw data with prediction data
// double final_pred_year = prediction year as inputed by user
// int month = prediction month as inputed by user
{
	// do the fitting
	string fitting_file_name = choice_two(data, final_pred_year);

	// load the fitting data into an array
	matrix pred;
	load_data(fitting_file_name, 1, pred);

	// workout the error at every point
	double const_error = error_overall(data, pred);

	// create the max and min arrays
	//now i want to get forecast_T +- error at each year:
	matrix forecast_upper = pred;
	matrix forecast_lower = pred;

	for (int i = 0; i < pred.n_row(); i++)
	{
		forecast_upper.set(i, 1, forecast_upper.i(i, 1) + const_error);
		forecast_lower.set(i, 1, forecast_lower.i(i, 1) - const_error);
	}
	string upper_file = save_to_csv(forecast_upper, "upper_");
	string lower_file = save_to_csv(forecast_lower, "lower_");

	string plot_command = "plot '" + fitting_file_name + "' using 1:2 with lines title 'Regression', '" +
						  original_file + "' every ::1 using 1:" +
						  to_string(month + 1) + " with points pointtype 6 title 'Data', '" + upper_file +
					      "' using 1:2 with lines title 'Upper Bound', '" + lower_file + "' using 1:2 with lines title 'Lower Bound'\n";
	plot(plot_command, "Tokyo Linear Regression with Time Invariant Error");

	// delete unwanted files
	remove(upper_file.c_str());
	remove(lower_file.c_str());
	remove(fitting_file_name.c_str());
	return;
}


// workout error increasing in time
void choice_five(matrix &data_trunc, matrix &data_orig,
				 double final_pred_year, int month, string original_file)
// matrix &data_trunc = temperature and years within model range (at specified month)
// matrix &data_orig = temperature and years of input file eg JMA_tokyo_data.csv
// double final_pred_year = prediction year as inputed by user
// int month = prediction month as inputed by user
// string original file = eg JMA_tokyo_data.csv, to plot raw data with prediction data
{
	// do the fitting
	string fitting_file_name = choice_two(data_trunc, final_pred_year);

	// load the fitting data into an array
	matrix data_pred(2), error(2);
	load_data(fitting_file_name, 1, data_pred);

	// compute last year of model from truncated year vector
	double model_data_end = 0;
	model_data_end = data_trunc.from_end(0, 0);
	error_per_point(data_pred, data_orig, error, model_data_end);

	// fitting to the errors and saving
	string error_fit_file_name = choice_two(error, final_pred_year, "d_error");
	
	// loading that save
	matrix pred_error(2);
	load_data(error_fit_file_name, 1, pred_error);

	// create the max and min arrays to get forecast_T +- error at each year:
	matrix forecast_upper = pred_error;
	matrix forecast_lower = pred_error;
	double diff = pred_error.i(0, 0) - data_pred.i(0,0);

	for (double i = 0; i < pred_error.n_row(); i++)
	{
		forecast_upper.set(i, 1, data_pred.i(i + diff, 1) + forecast_upper.i(i, 1));
		forecast_lower.set(i, 1, data_pred.i(i + diff, 1) - forecast_lower.i(i, 1));
	}

	string upper_file = save_to_csv(forecast_upper, "upper_");
	string lower_file = save_to_csv(forecast_lower, "lower_");

	string plot_command = "plot '" + upper_file + "' using 1:2 with lines title 'Upper Bound', '" +
						   lower_file + "' using 1:2 with lines title 'Lower Bound', '" +
						   fitting_file_name + "' using 1:2 with lines title 'Regression', '" +
						   original_file + "' every ::1 using 1:" + to_string(month + 1) + " with points pointtype 17 title 'Data'\n";

	plot(plot_command, "Tokyo Linear Regression with Time Dependent Error");

	return;
}


void menu()
{
	// input vairables
	string filepath;
	char repeat;
	int user_choice, pred_month;
	double start_year, end_year, pred_year;
	matrix trunc_data, orig_data;

	cout << "\n==========================\n\tMenu\n==========================\n";


	//user choice
	cout << "Command \t Description"
		<< "\n --------------------------------------------------------------------------"
		<< "\n 1 \t\t Temperature prediction for one specified year using "
		<< "\n    \t\t  both linear regression and lagrange interpolation"
		<< "\n 2 \t\t Prediction up to a specified year using linear regression."
		<< "\n     \t\t  Data will be saved locally in .csv file automatically"
		<< "\n 3 \t\t Prediction up to a specified year. Linear regression "
		<< "\n     \t\t  is used. Data is saved and plotted"
		<< "\n 4 \t\t Prediction up to a specified year with time invariant uncertanties. "
		<< "\n     \t\t  Linear regression is used. Data is plotted"
		<< "\n 5 \t\t Prediction up to a specified year with time varying uncertanties "
		<< "\n     \t\t  found with linear regression of error on unused data points."
		<< "\n     \t\t  Data is saved and plotted"
		<< "\n\nEnter a command: ";
	user_choice = input<int>(5, 1);


	// input file path
	cout << "\n\nDefault files for Tokyo climate data are: \n JMA_tokyo_data.csv\n"
		<< " JMA_tokyo_max.csv \n"
		<< " JMA_tokyo_min.csv"
		<< "\nPlease input the path to the climate data of your choice (.csv): ";
	cin >> filepath;
	//filepath = "JMA_tokyo_data.csv";


	// input truncation dates
	cout << "\nHistorical Data for Tokyo temperature goes back to 1876."
		<< "\nPlease write the year at which you would like to start the model: ";
	start_year = input<double>(2013, 1876);
	cout << "\nEnd year of model: ";
	end_year = input<double>(2018, start_year);


	// input prediction month
	map<const int, string> map_month;
	init_month(map_month);
	cout << "\n Please indicate the month you'd like to create the model for"
		<< "\n(JAN: 1, FEB: 2, MAR: 3, APR: 4, MAY: 5, JUN: 6, JUL: 7, AUG : 8,\n"
		<< "SEP : 9, OCT : 10, NOV : 11, DEC : 12, Year Average : 13)\n"
		<< "Prediction month: ";
	pred_month = input<int>(13, 1);
	cout << "\nSummary of model: \n" << "Model created from historical data between "
		<< start_year << " and " << end_year << " for " << map_month.find(pred_month)->second << "\n\n";


	// input year to predict up to if needed
	if (user_choice != 1)
	{
		cout << "Wanted projection final year (capped at 2100): ";
		pred_year = input<double>(2100, start_year);
	}

	// load and truncate the data
	load_data(filepath, pred_month, orig_data);
	truncate_data(trunc_data, start_year, end_year, orig_data);

	// select the correct program
	switch (user_choice)
	{
	case 1:
		choice_one(trunc_data, pred_month);
		break;
	case 2:
		choice_two(trunc_data, pred_year);
		break;
	case 3:
		choice_three(trunc_data, pred_month, filepath, pred_year);
		break;
	case 4:
		choice_four(trunc_data, filepath, pred_year, pred_month);
		break;
	case 5:
		choice_five(trunc_data, orig_data, pred_year, pred_month, filepath);
		break;
	}
}


int main()
{
	// input vairables
	string filepath;
	char repeat;
	int user_choice, pred_month;
	double start_year, end_year, pred_year;
	matrix trunc_data, orig_data;

	// intro
	cout << "Climate change predictor written by Deirdree Polak and Richard Boyne\n"
		 << "for ACSE 5.1 coursework - Imperial College London\n";
	
	while (true)
	{
		// call the menu every iteration
		menu();

		// restart program if wanted
		cout << "\n\nEnter q to quit or anyother key to restart: ";
		cin >> repeat;
		if (repeat == 'q') return 0;

	}
}

