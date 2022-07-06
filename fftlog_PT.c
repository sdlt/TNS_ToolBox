#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fftw3.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "fftlog_PT.h"

/*----------------------------------------------------------------*
 *FFTLog PT                                                       *
 *----------------------------------------------------------------*/

void FFTLogPT(FFTLog_config *fc, double (*func)(double))
{
  int i;
  double L=log(fc->max/fc->min),dlogk = log(fc->max/fc->min)/(double)fc->N;

  for(i=0;i<fc->N;i++){
    double k = fc->min*exp((double)i*dlogk);
    fc->an[i][0] = (*func)(k)*pow(k,-fc->q);
    fc->an[i][1] = 0;
  }

  fftw_execute(fc->p_forward);

  for(i=1;i<fc->N/2+1;i++){
    gsl_complex cm = gsl_complex_rect(fc->cm[i][0],fc->cm[i][1]);
    gsl_complex c = gsl_complex_exp(gsl_complex_mul_real(gsl_complex_rect(0,-2*M_PI*(double)i/L),log(fc->min)));
    gsl_complex m = gsl_complex_mul(cm,c);
    fc->cm[i][0] = GSL_REAL(m);
    fc->cm[i][1] = GSL_IMAG(m);
    fc->cm[fc->N-i][0] = fc->cm[i][0];
    fc->cm[fc->N-i][1] = -fc->cm[i][1];
  }
}

FFTLog_config *FFTLogPT_init(int N, double min, double max, double nu)
{ 
  FFTLog_config *fc = (FFTLog_config*)malloc(sizeof(FFTLog_config));
  double L;
  int i;

  fc->min        = min;
  fc->max        = max;
  fc->q          = nu;
  fc->mu         = 0;
  fc->N          = N;
  fc->an         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->cm         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->um         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  
  fc->p_forward = fftw_plan_dft_1d(N,fc->an,fc->cm,FFTW_FORWARD,FFTW_ESTIMATE);

  L = log(fc->max/fc->min);
  
  fc->um[0][0] = -0.5*fc->q;
  fc->um[0][1] = 0;
  for(i=1;i<fc->N/2+1;i++){
    fc->um[i][0] = -0.5*fc->q;
    fc->um[i][1] = -M_PI*(double)i/L;
    fc->um[fc->N-i][0] = fc->um[i][0];
    fc->um[fc->N-i][1] = -fc->um[i][1];
  }

  return fc;
}

void FFTLogPT_free(FFTLog_config *fc)
{  
  fftw_destroy_plan(fc->p_forward);
  fftw_free(fc->an);
  fftw_free(fc->cm);
  fftw_free(fc->um);  
  free(fc);

  return;
}
