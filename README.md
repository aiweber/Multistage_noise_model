# Multistage_noise_model

This repository contains Matlab code to fit a multiply stochastic model to neural data.

To get started, open `demo_fit_multistage_noise_model.m`.



In order to run this code, you will need:

1) The open-source Chebfun package, used for function approximation with Chebyshev polynomials.  It can be found at ([https://www.chebfun.org/](https://www.chebfun.org/))

2) A compiled mex file of ‘mult_likelihood.cpp’.  This requires the GNU Scientific Library (GSL), a free numerical library for C and C++.  The library and documentation can be found at ([https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)).


Although not necessary, this code is best run with MATLAB’s Parallel Computing Toolbox on a multi-core computer or computer cluster.

Additionally, this package includes the 'fminsearchbnd' function written by John D'Errico, which imposes bounds on MATLAB's built-in 'fminsearch' function.  Additional information can be found ([here](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon)).

