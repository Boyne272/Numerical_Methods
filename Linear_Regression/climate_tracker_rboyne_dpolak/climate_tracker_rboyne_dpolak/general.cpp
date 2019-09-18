#include "header.h"
#include "Matrix.h"


void load_data(string filepath, int column, matrix &mat)
//string filepath = filepath of original data file eg JMA_tokyo_data.csv
//int column = for  eg JMA_tokyo_data.csv, nth column is the nth month inputed by user
// matrix &mat = to store temperature and year data from data file eg JMA_tokyo_data.csv
{
	// check input matrix is empty
	assert(mat.size == 0);
	mat.n_col = 2;

	// create the file strem
	ifstream myFile;
	myFile.open(filepath, fstream::in);
	assert(!myFile.fail());

	// set up result stores
	string line, element;

	// skip the first line
	getline(myFile, line, '\n');

	// loop over lines
	while (getline(myFile, line, '\n'))
	{
		// create a stream for each line
		istringstream line_stream(line);

		// for every column
		for (int i = 0; i < 14; i++)
		{
			// get commer sperated string
			getline(line_stream, element, ',');
			// store if wanted
			if (i == 0)
				mat.append(stod(element));
			else if (i == column)
				mat.append(stod(element));
		}
	}

	// close the file
	myFile.close();
}


string get_current_time()
// get current time for file name
{
	char tmp_char[30];

	// method 1
	/*auto raw_time = chrono::system_clock::now();
	auto time_obj = chrono::system_clock::to_time_t(raw_time);
	ctime_s(tmp_char, sizeof tmp_char, &time_obj);*/

	// method 2
	time_t now = time(0);
	struct tm time_obj;
	localtime_s(&time_obj, &now);
	strftime(tmp_char, sizeof tmp_char, "%m-%d_%H%M%S", &time_obj);

	// return the strings
	string time_string = tmp_char;
	return time_string;
}


string save_to_csv(matrix &mat, string prefix = "data_", string header = "")
// matrix &mat = temperature and year data matrix
// string prefix = (optional) .csv file name prefix
// string header = (optional) header of .csv file columns
{
	// check both have same size
	assert(mat.filled);

	// get current time for file name
	string filename = prefix + get_current_time() + ".csv";

	//create output file
	ofstream save_file;
	save_file.open(filename);

	// write the header
	if (header == "")
	{
		for (int n = 0; n < mat.n_col; n++)
			save_file << "Col " << to_string(n) << ",";
		save_file << endl;
	}
	else
		save_file << header << endl;

	// loop over the rows and store the values within
	for (int i = 0; i < mat.n_row(); i++)
	{
		for (int j = 0; j < mat.n_col; j++)
			save_file << to_string(mat.i(i, j)) << ",";
		save_file << endl;
	}

	// close the file and exit
	save_file.close();
	cout << "\ndata saved to: " << filename;
	return filename;
}


void truncate_data(matrix &data_trunc, double &start_year, double &end_year,
				   matrix &data_orig)
// matrix &data = to store temperature and years within model range (at specified month)
// double &start_year, double &end_year = model range
// matrix &data_orig = temperature and years of input file eg JMA_tokyo_data.csv (at specified month)
{
	// check the years are in the bounds
	assert(start_year >= data_orig.i(0, 0));
	assert(end_year <= data_orig.from_end(0, 0));

	// check every year has a temperature
	assert(data_orig.filled);

	// check result matrix is empty
	assert(data_trunc.size == 0);
	data_trunc.n_col = 2;

	for (int i = 0; i < data_orig.n_row(); i++)
	{
		// if between the two years wanted
		if (end_year >= data_orig.i(i, 0) && data_orig.i(i, 0) >= start_year)
		{
			data_trunc.append(data_orig.i(i, 0));
			data_trunc.append(data_orig.i(i, 1));
		}
	}

	// check the matrix is correctly filled
	assert(data_trunc.filled);
}

void plot(string plot_string, string title = "Tokyo",
	string xlabel = "Time (years)", string ylabel = "Temp (Deg.C)")
// string plot_string = command line for GNU plot, passed in
// string title, string xlabel, string ylabel = plot title and axis names
{
	// create gnuplot command pipe and check it exists
	FILE *gpPipe = _popen("gnuplot -persist", "w");
	assert(gpPipe);

	string tmp_string;

	// setup headings
	tmp_string = "set title '" + title + "'\n";
	fprintf(gpPipe, tmp_string.c_str());
	tmp_string = "set xlabel '" + xlabel + "'\n";
	fprintf(gpPipe, tmp_string.c_str());
	tmp_string = "set ylabel '" + ylabel + "'\n";
	fprintf(gpPipe, tmp_string.c_str());

	// other settings
	fprintf(gpPipe, "set datafile separator ','\n");
	fprintf(gpPipe, "set pointsize 1\n");
	fprintf(gpPipe, "set key top left\n");
	//fprintf(gpPipe, "set terminal wxt size 500,500");


	// plot what is wanted
	fprintf(gpPipe, plot_string.c_str());

	// run the code in the pipe; exit and close pipe
	fflush(gpPipe);
	fprintf(gpPipe, "\nexit \n");
	_pclose(gpPipe);
}


void init_month(map<const int, string> &map_month)
// map<const int, string> &map_month = empty map to initialize
{
	map_month.insert(std::pair<const int, std::string>(1, "January"));
	map_month.insert(std::pair<const int, std::string>(2, "February"));
	map_month.insert(std::pair<const int, std::string>(3, "March"));
	map_month.insert(std::pair<const int, std::string>(4, "April"));
	map_month.insert(std::pair<const int, std::string>(5, "May"));
	map_month.insert(std::pair<const int, std::string>(6, "June"));
	map_month.insert(std::pair<const int, std::string>(7, "July"));
	map_month.insert(std::pair<const int, std::string>(8, "August"));
	map_month.insert(std::pair<const int, std::string>(9, "September"));
	map_month.insert(std::pair<const int, std::string>(10, "October"));
	map_month.insert(std::pair<const int, std::string>(11, "November"));
	map_month.insert(std::pair<const int, std::string>(12, "December"));
	map_month.insert(std::pair<const int, std::string>(13, "Average Year"));
}


