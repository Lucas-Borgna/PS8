                                                         Glen Cowan
                                                         RHUL Physics
                                                         27 November, 2013

                A sample fitting macro for root

simpleFit.C is a simple root macro to read in data from a file, do a
least-squares fit and make a plot.  It can be used as a starting point
for creating further fit macros.  To run, first run root and then

.L simpleFit.C
simpleFit()

and then enter the data file name when prompted.

Data are read in from a file (testData.txt) with 3 columns of numbers
representing x, y and dy.  Lines can be commented out with a #, and
the routine can easily be generalized to deal with a different number
of columns.

The fit function is contained at the top of the file simpleFit.C, and
is called fitFunc.  Currently this is a polynomial of order npar-1
(i.e., npar free parameters).  The value of npar is set as a global
variable at the top of the macro.

The output is stored in a gif and also eps file.  The macro allows
user access to the fit results (parameter values, errors, covariance
matrix, chi2).  The covariance matrix is stored in an object of type
TMatrixDSym, and this class is used also to find its inverse.



