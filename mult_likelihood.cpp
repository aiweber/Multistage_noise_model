#define _USE_MATH_DEFINES /* includes pi */
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "gsl/gsl_integration.h"
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "mex.h"
#include <cmath>


/* Computes pdf of responses by taking integral over multiplicative noise. */



double pfunc_mult (double k, void *params)
{
  double *alpha = (double *)params;
  double gObs = *(alpha);
  double upNoise = *(alpha+1);
  double multNoise = *(alpha+2);
  double y = *(alpha+3);
  double a = *(alpha+4);
  double b = *(alpha+5);
  double c = *(alpha+6);
  double d = *(alpha+7);
  
  
  double expCut = 500;
  double x = 0;
  double fl = 0;
  
  if((k-d)/a >= expCut){
     x = ((k-d)/a-c)/b;
     fl =  1 / (upNoise * sqrt(2*M_PI)) * exp(-pow(x-gObs,2)/(2*pow(upNoise,2))) / (a*b);
  }
  else{
     x = (log(exp((k-d)/a)-1)-c)/b;
     fl =  1 / (upNoise * sqrt(2*M_PI)) * exp(-pow(x-gObs,2)/(2*pow(upNoise,2))) * exp((k-d)/a) / (a*b*(exp((k-d)/a)-1));
  }

return fl / (multNoise*sqrt(2*M_PI*k)) * exp(-pow(y-k,2)/(2*k*pow(multNoise,2)));

}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
     #define y_out plhs[0]
     double *params = mxGetPr(prhs[0]);
     double res, err;

     double *xl = mxGetPr(prhs[1]);
     double *xu = mxGetPr(prhs[2]);

    
    // QAG integration
    
     int key = 6;        // indicates number of points on each integration step

     gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    
     gsl_function F;
     F.function = &pfunc_mult;
     F.params = params;
   
    gsl_integration_qag (&F, *(xl), *(xu), 1e-4, 1e-4, 1000, key, w, &res, &err);
    gsl_integration_workspace_free (w);

    
    y_out = mxCreateDoubleScalar(res);

}


