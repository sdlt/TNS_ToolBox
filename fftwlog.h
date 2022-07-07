#include <fftw3.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>

typedef struct NcmFftlog
{
  int Nr;
  int N;
  int N_2;
  int Nf;
  int Nf_2;
  size_t nderivs;
  size_t pad;
  double lnk0;
  double lnr0;
  double Lk;
  double Lk_N;
  double pad_p;
  double q;
  double ell;
  gsl_vector *lnr_vec;
  void **Gr_vec;
  fftw_complex *Fk;
  fftw_complex *Cm;
  fftw_complex *Gr;
  fftw_complex *CmYm;
  void **Ym;
  fftw_plan p_Fk2Cm;
  fftw_plan p_CmYm2Gr;
} NcmFftlog;

gsl_vector *vector_new (size_t n);
int fftwlog_get_full_size (NcmFftlog *fftlog);
double fftwlog_get_length (NcmFftlog *fftlog);
double fftwlog_get_full_length (NcmFftlog *fftlog);
double fftwlog_get_ell (NcmFftlog *fftlog);
double fftwlog_get_lnk0 (NcmFftlog *fftlog);
double fftwlog_get_lnr0 (NcmFftlog *fftlog);
int fftwlog_get_mode_index (NcmFftlog *fftlog, int i);
int fftwlog_get_array_index (NcmFftlog *fftlog, int phys_i);
void fftwlog_set_best_lnr0 (NcmFftlog *fftlog);
void fftwlog_set_best_lnk0 (NcmFftlog *fftlog);
NcmFftlog *fftwlog_init (int nsample, double min, double max, double q, double ell);
size_t util_fact_size (size_t n);
void fftwlog_free_all (NcmFftlog *fftlog);
void fftwlog_set_size (NcmFftlog *fftlog, int n);
void fftwlog_get_Ym(NcmFftlog *fftlog, fftw_complex *Ym_0);
void fftwlog_eval_by_gsl_function (NcmFftlog *fftlog, gsl_function *Fk);
void fftwlog_eval (NcmFftlog *fftlog);
gsl_vector *fftwlog_output (NcmFftlog *fftlog, int id);
int xi_fftw_vec(gsl_function *Pk, double *vec, int nr, double rlog0, double drlog, NcmFftlog *fftlog, int mode);
int xi_fftw_vec_der(gsl_function *Pk, double *vec, double *dvec, int nr, double rlog0, double drlog, double alpha, NcmFftlog *fftlog);
int xi_fftw_vec_lin(gsl_function *Pk, double *vec, int nr, double rmin, double rmax, NcmFftlog *fftlog, int mode);
int xi_fftw_vec_der_lin(gsl_function *Pk, double *vec, double *dvec, int nr, double rmin, double rmax, double alpha, NcmFftlog *fftlog);
