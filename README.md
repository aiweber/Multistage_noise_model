# Multistage_noise_model

Associated with ([Weber, Shea-Brown, & Rieke
2020](https://www.biorxiv.org/content/10.1101/2020.02.17.951830v1)), this repository contains MATLAB code to fit a multiply stochastic model to neural data.

To get started, open `demo_fit_multistage_noise_model.m`.



In order to run this code, you will need:

1) The open-source Chebfun package, used for function approximation with Chebyshev polynomials.  It can be found at ([https://www.chebfun.org/](https://www.chebfun.org/))

2) A compiled mex file of `mult_likelihood.cpp`.  This requires the GNU Scientific Library (GSL), a free numerical library for C and C++.  The library and documentation can be found at ([https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)).


Although not necessary, this code is best run with MATLAB’s Parallel Computing Toolbox on a multi-core computer or computer cluster.

Additionally, this package includes the `fminsearchbnd` function written by John D'Errico, which imposes bounds on MATLAB's built-in fminsearch function.  Additional information can be found ([here](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon)).



A note to Windows users on compiling the mex file: 
The mex file can be more difficult to compile on a Windows machine than on Mac or Linux.  Keep in mind that you should use the same compiler to compile both the GNU Scientific Library and the mex file.  The author has had success using MinGW-MSYS (not MSYS2).  First, install MinGW.  Then, use the mingw-get package to install MSYS: `mingw-get install msys`.  Within MATLAB, you may also need to run the `configuremingw` script (available in MATLAB Central) to get MATLAB to recognize your installation of MinGW.