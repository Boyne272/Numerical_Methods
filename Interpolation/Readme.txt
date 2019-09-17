Richard Boyne
CID 01057503
Submitted 13/11/17 at 11am with 

------------------------------------------------

Task 1 - folating_point_vairables.py:

Run program then call solve_f() 

no float-128 as university computers automatically case these as float-64

n.b. solve is also the validate function here

------------------------------------------------

Task 2 - matrix_methods.py:

Run program then call solve_LU() to find solutions and automatically check them

call validate_LU() to run for a 100 by 100 matrix with values from 0 to 10 and
then check the solution found (can change the size, seed and range for the matrix
with parameters size, r, s).

------------------------------------------------

Task 3 - interpoilation:

Run program then call solve_i() to find solutions

call validate_i() to run for 20 y random points with uniform x spacing or with
xrand = True for random c and y points (other options are the number, seed and
range of points with n,s,r respectivly)

------------------------------------------------

Task 4 - fourier_transforms.py:

Just run program for three plots:

convolution both from numpy for reference and via the inverse FFT product
the original functions and their FFTs for reference
the inverse FFTs and fourier transform products for reference
the last two can be turned off by hashing line 95 (the final line)


Supplement - extra_fft_attempt.py:

Although this was not requiered it help understanding and so I thought I would
include it. it is just a self written fft that gives the same output as the numpy fft for
a box function

Run the program to see the comparison with a dft and numpy fft
Note this section is undercommented as I did not have time and it is not requiered
code

------------------------------------------------

Task 5 - random_numbers.py:

Just run program for three plots that are the three random number distrubutions

The verify function is automatically called and prints the mean square difference
of each point to the expected pdfs (goes to zero if n goes to inf)

the number of bins and number of values generated can be changed with res and n
respectivly on lines 69,68