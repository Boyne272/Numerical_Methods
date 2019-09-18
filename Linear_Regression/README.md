# acse-5-advanced-programming-coursework-1-tokyo

By Richard Boyne and Deirdree Polak for ACSE 5 coursework 1

VS project file is found in "acse-5-advanced-programming-coursework-1-tokyo/climate_tracker_rboyne_dpolak/" with data and gnuplot.exe in "acse-5-advanced-programming-coursework-1-tokyo/climate_tracker_rboyne_dpolak/climate_tracker_rboyne_dpolak/" (i.e. the local directory for code run in VS strudio IDE). Standalone executables can be found in dubug if requiered.

## Code Requierements
Pre-requisites detailed in the "header.h" file are all part of the C++ STL so do not require additional installation. For any of the plotting options (i.e. opions 3, 4 and 5) gnuplot must be installed with all "*.dll" files in the computers "/system_32/" directory and gnuplot.exe must be in the same directory as the executable file (included in github repository). If you don't wish to install these see option 2 below.

## Data input
The code inputs a "*.csv" file with the first column being year of a measurement and then multiple columns for different months (with an additional column for the annual average if wanted). In this repository are 3 csv files sourced from Japan Meteorological Agency (1-3-4 Otemachi, Chiyoda-ku, Tokyo 100-8122, Japan, https://www.data.jma.go.jp/obd/stats/data/en/smp/index.html) giving complete records for temperature of Tokyo city recorded to 0.1 degree celsius dating back to 1876. These are stored in the following files
  
  - JMA_tokyo_data.csv has the monthly average (and annual) temperature
  - JMA_tokyo_max.csv has the monthly maximum (and annual) temperature 
  - JMA_tokyo_min.csv has the monthly minimum (and annual) temperature
  
## User Interface
Upon execution the user is prompted with a menu where they may select a desired program from the following options:
 1               Temperature prediction for one specified year using
                  both linear regression and lagrange interpolation
 2               Prediction up to a specified year using linear regression.
                  Data will be saved locally in .csv file automatically
 3               Prediction up to a specified year. Linear regression
                  is used. Data is saved and plotted
 4               Prediction up to a specified year with time invariant uncertanties.
                  Linear regression is used. Data is plotted
 5               Prediction up to a specified year with time varying uncertanties
                  found with linear regression of error on unused data points.
                  Data is saved and plotted

After selecting the user will be prompted for:
- a csv file to source data from
- a start and end year for the model to be based on
- a month of the year to predict with (choosing month 13 gives the annual average)
- a year to predict for/up to. 

After this the individual program will run, and once done the option is given to return to the menu or to exit the program.

## Data Output
Program 1 saves no data and simply prints to the user in the terminal.

Programs 2, 3 and 5 save data as .csv files in the same directory as the executable. These files are given a prefix dependent on their type:
- data_ stores regression predictions
- upper_ and lower_ stores error bound estimations
- derror stores time dependent error estimations
each prefix is then followed by the date and time of creation to prevent overwriting of pre-existing results.

Programs 3, 4 and 5 also plot results in gnuplot, these can be saved directly from the gnuplot GUI that appears on execution.

## Functions & Classes
### class Matrix - Matrix.h

The Matrix class is essential to the running of the program. It creates a data structure for the loaded data that is then replicated for the processed data. It consists of 9 methods, including matrix multiplication, transpose and inversion. It also deals with operations such as append and looking up the last value.

### load_data and truncatedata -- general.cpp

The data is loaded using fstream. The data is loaded up from the file path inputted by user. The data corresponding to the month of prediction is saved in a matrix. The data is then truncated to fit wanted model.

### lagrange and linearregress -- numerical.cpp
The temperature prediction is computed through lagrange and linear regression for choice one, and then only by linear regression for choice two to five. Standard deviation is added to choice four, and the increasing error is computed for choice five.

### save_to_csv and get_time -- general.cpp
This opens a fstream to save the temperature data on a .csv file. The other functions is used to have a timestamp on the name of the .csv file. 

### plot -- general.cpp
GNU Plot function. It is called by choice 3, 4 and 5. The function pipes to the gnu command line to create a plot that can be saved by the user.

### input -- climate_tracker_rboyne_dpolak.cpp
User Input function to limit the type and value range of user input.

### menu -- climate_tracker_rboyne_dpolak.cpp
Function to bring up user menu. This allows for the restart functionality in the main() to go back to the start menu or quit program.
