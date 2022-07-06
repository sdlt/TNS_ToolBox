#include "fftwlog.h"
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <omp.h>

#define fftw_alloc_complex(n) (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * (n))

gsl_vector *vector_new (size_t n)
{
  gsl_vector *d = gsl_vector_alloc(sizeof(double)*n);
  return d;
}

int fftwlog_get_full_size (NcmFftlog *fftlog)
{
  return fftlog->Nf;
}

double fftwlog_get_length (NcmFftlog *fftlog)
{
  return fftlog->Lk;
}

double fftwlog_get_full_length (NcmFftlog *fftlog)
{
  return fftlog->Lk + 2.0 * fftlog->Lk_N * fftlog->pad;
}

double fftwlog_get_ell (NcmFftlog *fftlog)
{
  return fftlog->ell;
}

double fftwlog_get_lnk0 (NcmFftlog *fftlog)
{
  return fftlog->lnk0;
}

double fftwlog_get_lnr0 (NcmFftlog *fftlog)
{
  return fftlog->lnr0;
}

int fftwlog_get_mode_index (NcmFftlog *fftlog, int i)
{
  return (i > fftlog->Nf_2) ? i - fftlog->Nf : i;
}

int fftwlog_get_array_index (NcmFftlog *fftlog, int phys_i)
{
  return (phys_i < 0) ? phys_i + fftlog->Nf : phys_i;
}

void fftwlog_set_best_lnr0 (NcmFftlog *fftlog)
{
  double lnk0 = fftwlog_get_lnk0 (fftlog);
  double Lk   = fftwlog_get_length (fftlog);
  double ell  = fftwlog_get_ell (fftlog);
  double lnc0 = (fabs(ell)<=1e-12) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - log(M_PI) + 2.0 * gsl_sf_lngamma (1.5 + ell)) / (2.0 * (1.0 + ell));
  
  fftlog->lnr0 = -lnk0; //+lnc0;
}

void fftwlog_set_best_lnk0 (NcmFftlog *fftlog)
{
  double lnr0 = fftwlog_get_lnr0 (fftlog);
  double Lk   = fftwlog_get_length (fftlog);
  double ell  = fftwlog_get_ell (fftlog);
  double lnc0 = (fabs(ell)<=1e-12) ? 0.0 : ((ell - 1.0) * Lk + 2.0 * (ell + 1.0) * M_LN2 - log(M_PI) + 2.0 * gsl_sf_lngamma (1.5 + ell)) / (2.0 * (1.0 + ell));
  
  fftlog->lnk0 = -lnr0; //+lnc0;
}

NcmFftlog *fftwlog_init (int nsample, double min, double max, double q, double ell)
{
  NcmFftlog *fftlog = (NcmFftlog*)malloc(sizeof(NcmFftlog));
  
  fftlog->lnr0  = 0.0;
  fftlog->lnk0  = log(max*min)*0.5;
  fftlog->Lk    = log(max/min);
  fftlog->Lk_N  = 0.0;
  fftlog->pad_p = 0.2;
  fftlog->Nr    = 0;
  fftlog->N     = 0;
  fftlog->N_2   = 0;
  fftlog->Nf    = 0;
  fftlog->Nf_2  = 0;
  fftlog->pad   = 0;

  fftlog->nderivs = 0;
  fftlog->q       = q;
  fftlog->ell     = ell;

  fftlog->lnr_vec   = NULL;
  fftlog->Gr_vec    = NULL;
  fftlog->Ym        = NULL;
  
  fftlog->Fk        = NULL;
  fftlog->Cm        = NULL;
  fftlog->Gr        = NULL;
  fftlog->CmYm      = NULL;
  fftlog->p_Fk2Cm   = NULL;
  fftlog->p_CmYm2Gr = NULL;
  
  fftlog->Gr_vec = malloc(sizeof(void*)*(fftlog->nderivs+1));
  fftlog->Ym  = malloc(sizeof(void*)*(fftlog->nderivs+1));
  
  fftwlog_set_best_lnr0 (fftlog);
  fftwlog_set_size (fftlog, nsample);
  
  return fftlog;
}
 
size_t util_fact_size (size_t n)
{
  if (n == 1) return 0;
  else {
   size_t r2 = n % 2;
   size_t r3 = n % 3;
   size_t r5 = n % 5;
   size_t r7 = n % 7;
   size_t m = 1;
  
   if (r2 == 0) m *= 2;
   if (r3 == 0) m *= 3;
   if (r5 == 0) m *= 5;
   if (r7 == 0) m *= 7;
  
   if (m != 1) {
    if (n / m == 1) return m;
    else return m * util_fact_size (n / m);
   } else return util_fact_size (n + 1);
  }
}

void fftwlog_free_all (NcmFftlog *fftlog)
{
  int i;
  
  fftw_free(fftlog->Fk);
  fftw_free(fftlog->Cm);
  fftw_free(fftlog->CmYm);
  fftw_free(fftlog->Gr);

  fftw_destroy_plan(fftlog->p_Fk2Cm);
  fftw_destroy_plan(fftlog->p_CmYm2Gr);

  gsl_vector_free(fftlog->lnr_vec);

  for (i=0; i<fftlog->nderivs+1; i++) {
    gsl_vector_free(fftlog->Gr_vec[i]);
    fftw_free(fftlog->Ym[i]);
  } 
  free(fftlog->Gr_vec);
  free(fftlog->Ym);
}
 
void fftwlog_set_size (NcmFftlog *fftlog, int n)
{
  size_t nt = n * (1.0 + fftlog->pad_p);

  fftlog->Nr = n;

  nt = util_fact_size (nt);
  fftlog->pad = nt * fftlog->pad_p * 0.5 / (1.0 + fftlog->pad_p);
  n = nt - 2 * fftlog->pad;
  
  if ((n != fftlog->N) || (n + 2 * fftlog->pad != fftlog->Nf))
  {
    size_t i;
    
    fftlog->N    = n;
    fftlog->N_2  = fftlog->N / 2;
    fftlog->Lk_N = fftlog->Lk / (1.0 * fftlog->N);

    fftlog->Nf   = fftlog->N + 2 * fftlog->pad;
    fftlog->Nf_2 = fftlog->N_2 + fftlog->pad;

    fftlog->Fk      = fftw_alloc_complex (fftlog->Nf);
    fftlog->Cm      = fftw_alloc_complex (fftlog->Nf);
    fftlog->CmYm    = fftw_alloc_complex (fftlog->Nf);
    fftlog->Gr      = fftw_alloc_complex (fftlog->Nf);

    fftlog->lnr_vec = vector_new (fftlog->N);

    fftlog->p_Fk2Cm   = fftw_plan_dft_1d (fftlog->Nf, fftlog->Fk, fftlog->Cm, FFTW_FORWARD, FFTW_DESTROY_INPUT);
    fftlog->p_CmYm2Gr = fftw_plan_dft_1d (fftlog->Nf, fftlog->CmYm, fftlog->Gr, FFTW_FORWARD, FFTW_DESTROY_INPUT);

    for (i=0; i<fftlog->nderivs+1; i++)
    {
      gsl_vector *Gr_vec_i = vector_new (fftlog->N);
      fftw_complex *Ym_i  = fftw_alloc_complex (fftlog->Nf);
      fftlog->Gr_vec[i] = Gr_vec_i;
      fftlog->Ym[i] = Ym_i;
    }
  }
}

void fftwlog_get_Ym(NcmFftlog *fftlog, fftw_complex *Ym_0)
{ 
  double twopi_Lt = 2.0 * M_PI / fftwlog_get_full_length (fftlog);
  int Nf          = fftwlog_get_full_size (fftlog);

  fftw_complex *Ym_base = (fftw_complex *) Ym_0;
  int i;
  
  double q = fftlog->q;
  double ell = fftlog->ell;

  if (q == 0.5)
  {
    for (i = 0; i < Nf; i++)
    {
      int phys_i = fftwlog_get_mode_index (fftlog, i);
      gsl_complex z, A, xup, two_x, U;
      gsl_sf_result lngamma_rho_up, lngamma_theta_up;

      A = gsl_complex_rect(0.5, twopi_Lt * phys_i);
      xup = gsl_complex_rect(0.5*(1.0 + ell + GSL_REAL(A)),0.5*GSL_IMAG(A));
      gsl_sf_lngamma_complex_e (GSL_REAL(xup), GSL_IMAG(xup), &lngamma_rho_up, &lngamma_theta_up);

      U = gsl_complex_polar(1,2*lngamma_theta_up.val);
      two_x = gsl_complex_pow(gsl_complex_rect(2,0), A);
      z = gsl_complex_mul(two_x, U);
      
      Ym_base[i][0] = GSL_REAL(z);
      Ym_base[i][1] = GSL_IMAG(z);
    }
  }
  else
  {
    for (i = 0; i < Nf; i++)
    {
      int phys_i = fftwlog_get_mode_index (fftlog, i);
      gsl_complex z, A, xup, xdw, two_x, U;
      gsl_sf_result lngamma_rho_up, lngamma_theta_up;
      gsl_sf_result lngamma_rho_dw, lngamma_theta_dw;

      A = gsl_complex_rect(q, twopi_Lt * phys_i);
      xup = gsl_complex_rect(0.5*(1.0 + ell + GSL_REAL(A)), 0.5*GSL_IMAG(A));
      xdw = gsl_complex_rect(0.5*(1.0 + ell - GSL_REAL(A)),-0.5*GSL_IMAG(A));
      gsl_sf_lngamma_complex_e (GSL_REAL(xup), GSL_IMAG(xup), &lngamma_rho_up, &lngamma_theta_up);
      gsl_sf_lngamma_complex_e (GSL_REAL(xdw), GSL_IMAG(xdw), &lngamma_rho_dw, &lngamma_theta_dw);
      
      U = gsl_complex_polar(exp(lngamma_rho_up.val-lngamma_rho_dw.val),lngamma_theta_up.val-lngamma_theta_dw.val); 
      two_x = gsl_complex_pow(gsl_complex_rect(2,0), A);
      z = gsl_complex_mul(two_x, U);
      
      Ym_base[i][0] = GSL_REAL(z);
      Ym_base[i][1] = GSL_IMAG(z);
    }
  }
}

void fftwlog_eval_by_gsl_function (NcmFftlog *fftlog, gsl_function *Fk)
{
  int i,dump;

  #pragma omp parallel for 
  for (i=0; i<fftlog->N; i++)
  {
    int phys_i   = i - fftlog->N_2;
    double lnk_i = fftlog->lnk0 + fftlog->Lk_N * phys_i;
    double k_i   = exp(lnk_i);
    double Fk_i  = Fk->function(k_i,Fk->params);

    fftlog->Fk[fftlog->pad + i][0] = Fk_i * pow(k_i,fftlog->q);
    fftlog->Fk[fftlog->pad + i][1] = 0;
  }
  
  fftwlog_eval (fftlog);
}

void fftwlog_eval (NcmFftlog *fftlog)
{
  size_t nd;
  int i;

  fftw_execute (fftlog->p_Fk2Cm);
  
  {
    double Lt          = fftwlog_get_full_length (fftlog);
    double twopi_Lt    = 2.0 * M_PI / Lt;
    fftw_complex *Ym_0 = fftlog->Ym[0];
    double lnr0k0      = fftlog->lnk0 + fftlog->lnr0;
    
    fftwlog_get_Ym (fftlog, Ym_0);

    for (i=0; i<5; i++) {
      double theta  = gsl_complex_arg(gsl_complex_rect(Ym_0[fftlog->Nf/2][0],Ym_0[fftlog->Nf/2][1]));
      double M      = (fftlog->Nf / Lt) * lnr0k0 - theta / M_PI;
      double long M_round  = M;
      double dM     = M - M_round;
      lnr0k0       -= (Lt / fftlog->Nf) * dM;
      fftlog->lnr0 -= (Lt / fftlog->Nf) * dM;
    }
    
    for (i=0; i<fftlog->Nf; i++)
    {
      int phys_i = fftwlog_get_mode_index (fftlog, i);
      gsl_complex z,arg;
      fftw_complex Ym_ndm1;

      arg = gsl_complex_polar (1, -lnr0k0 * twopi_Lt * phys_i); 
      z = gsl_complex_mul (gsl_complex_rect(Ym_0[i][0],Ym_0[i][1]), arg);
      Ym_0[i][0] = GSL_REAL(z);
      Ym_0[i][1] = GSL_IMAG(z);

      Ym_ndm1[0] = Ym_0[i][0];
      Ym_ndm1[1] = Ym_0[i][1];
      for (nd=1; nd<=fftlog->nderivs; nd++)
      {
        fftw_complex *Ym_nd = fftlog->Ym[nd];
	z = gsl_complex_rect(-Ym_ndm1[0], twopi_Lt * phys_i * Ym_ndm1[1]);
	Ym_ndm1[0] = Ym_nd[i][0] = GSL_REAL(z);
	Ym_ndm1[1] = Ym_nd[i][1] = GSL_IMAG(z);
      }
    }

    if ((fftlog->Nf % 2) == 0)
    {
      int Nf_2_index = fftwlog_get_array_index (fftlog, fftlog->Nf/2);
      
      for (nd=0; nd<=fftlog->nderivs; nd++)
      {
        fftw_complex *Ym_nd = fftlog->Ym[nd];
        Ym_nd[Nf_2_index][1] = 0;
      }
    }
  }

  for (i=0; i<fftlog->N; i++)
  {
    int phys_i = i - fftlog->N_2;
    double lnr = fftlog->lnr0 + phys_i * fftlog->Lk_N;
    gsl_vector_set (fftlog->lnr_vec, i, lnr);
  }
  
  for (nd=0; nd<=fftlog->nderivs; nd++)
  {
    double norm = (double)fftlog->Nf;
    gsl_vector *Gr_nd   = fftlog->Gr_vec[nd];
    fftw_complex *Ym_nd = fftlog->Ym[nd];
    
    for (i=0; i<fftlog->Nf; i++)
    {
      gsl_complex z = gsl_complex_mul(gsl_complex_rect(fftlog->Cm[i][0],fftlog->Cm[i][1]), gsl_complex_rect(Ym_nd[i][0],Ym_nd[i][1]));
      fftlog->CmYm[i][0] = GSL_REAL(z);
      fftlog->CmYm[i][1] = GSL_IMAG(z);
    }

    fftlog->CmYm[fftlog->Nf_2][1] = fftlog->CmYm[fftlog->Nf_2 + 1][1] = 0;

    fftw_execute (fftlog->p_CmYm2Gr);

    for (i=0; i<fftlog->N; i++)
    {
      double lnr_i   = gsl_vector_get (fftlog->lnr_vec, i);
      double r_i     = exp(lnr_i);
      double Gr_nd_i = fftlog->Gr[i + fftlog->pad][0] * pow(r_i,fftlog->q) / norm * pow(2*M_PI*r_i,-1.5);

      gsl_vector_set (Gr_nd, i, Gr_nd_i);
    }
  }
}

gsl_vector *fftwlog_output (NcmFftlog *fftlog, int id)
{
  gsl_vector *out;
  
  if (id<0 || id>fftlog->nderivs) id=0;
  out = (gsl_vector *)fftlog->Gr_vec[id];
  return out;
}

int xi_fftw_vec(gsl_function *Pk, double *vec, int nr, double rlog0, double drlog, NcmFftlog *fftlog, int mode)
{
  int i;
  double *logr = malloc(fftlog->N*sizeof(double));
  double *xir = malloc(fftlog->N*sizeof(double));
  gsl_vector *xi;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fftlog->N);
    
  if (mode<0 || mode>1) mode=0;

  fftwlog_eval_by_gsl_function (fftlog,Pk);
  xi = (gsl_vector *)fftlog->Gr_vec[0];
  
  for(i=0;i<fftlog->N;i++) {
    logr[i] = gsl_vector_get(fftlog->lnr_vec,i);
    xir[i]  = gsl_vector_get(xi,i);
  }
  gsl_spline_init(spline, logr, xir, fftlog->N);

  for (i=1;i<=nr;i++) {
    double rlog=rlog0+i*drlog;
    double r_prime=pow(10,rlog);
    if (!mode) vec[i]=gsl_spline_eval (spline, log(r_prime), acc);
    else vec[i]=gsl_spline_eval_deriv (spline, log(r_prime), acc);
  }

  free (logr);
  free (xir);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}

int xi_fftw_vec_der(gsl_function *Pk, double *vec, double *dvec, int nr, double rlog0, double drlog, double alpha, NcmFftlog *fftlog)
{
  int i;
  double *logr = malloc(fftlog->N*sizeof(double));
  double *xir = malloc(fftlog->N*sizeof(double));
  gsl_vector *xi;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fftlog->N);

  fftwlog_eval_by_gsl_function (fftlog,Pk);
  xi = (gsl_vector *)fftlog->Gr_vec[0];
  
  for(i=0;i<fftlog->N;i++) {
    logr[i] = gsl_vector_get(fftlog->lnr_vec,i);
    xir[i]  = gsl_vector_get(xi,i);
  }
  gsl_spline_init(spline, logr, xir, fftlog->N);

  for (i=1;i<=nr;i++) {
    double rlog=rlog0+i*drlog;
    double r_prime=pow(10,rlog);
    vec[i]=gsl_spline_eval (spline, log(r_prime*alpha), acc);
    dvec[i]=gsl_spline_eval_deriv (spline, log(r_prime*alpha), acc);
  }

  free (logr);
  free (xir);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}

int xi_fftw_vec_lin(gsl_function *Pk, double *vec, int nr, double rmin, double rmax, NcmFftlog *fftlog, int mode)
{
  int i;
  double *logr = malloc(fftlog->N*sizeof(double));
  double *xir = malloc(fftlog->N*sizeof(double));
  gsl_vector *xi;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fftlog->N);
  double dr = (rmax-rmin)/(double)nr;
    
  if (mode<0 || mode>1) mode=0;

  fftwlog_eval_by_gsl_function (fftlog,Pk);
  xi = (gsl_vector *)fftlog->Gr_vec[0];

  #pragma omp parallel for 
  for(i=0;i<fftlog->N;i++) {
    logr[i] = gsl_vector_get(fftlog->lnr_vec,i);
    xir[i]  = gsl_vector_get(xi,i);
  }
  gsl_spline_init(spline, logr, xir, fftlog->N);
  
  #pragma omp parallel for 
  for (i=1;i<=nr;i++) {
    double r_prime = rmin+(i-1)*dr+0.5*dr;
    if (!mode) vec[i]=gsl_spline_eval (spline, log(r_prime), acc);
    else vec[i]=gsl_spline_eval_deriv (spline, log(r_prime), acc);
  }

  free (logr);
  free (xir);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}

int xi_fftw_vec_der_lin(gsl_function *Pk, double *vec, double *dvec, int nr, double rmin, double rmax, double alpha, NcmFftlog *fftlog)
{
  int i;
  double *logr = malloc(fftlog->N*sizeof(double));
  double *xir = malloc(fftlog->N*sizeof(double));
  gsl_vector *xi;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fftlog->N);
  double dr = (rmax-rmin)/(double)nr;

  fftwlog_eval_by_gsl_function (fftlog,Pk);
  xi = (gsl_vector *)fftlog->Gr_vec[0];
  
  for(i=0;i<fftlog->N;i++) {
    logr[i] = gsl_vector_get(fftlog->lnr_vec,i);
    xir[i]  = gsl_vector_get(xi,i);
  }
  gsl_spline_init(spline, logr, xir, fftlog->N);

  for (i=1;i<=nr;i++) {
    double r_prime = rmin+(i-1)*dr+0.5*dr;
    vec[i]=gsl_spline_eval (spline, log(r_prime*alpha), acc);
    dvec[i]=gsl_spline_eval_deriv (spline, log(r_prime*alpha), acc);
  }

  free (logr);
  free (xir);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}
