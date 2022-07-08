#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <fftw3.h>
#include "fftlog_PT.h"
#include <omp.h>

FFTLog_config *fcpt;
gsl_matrix_complex *I_store;
double *kkk;

static gsl_interp_accel *acc[9];
static gsl_spline *spline[9];

extern double Pk_lin(double k);

gsl_complex Ifunc(int m1, int m2)
{
  gsl_sf_result g1,g2,g3,g4,g5,g6;
  gsl_sf_result ag1,ag2,ag3,ag4,ag5,ag6;
  double zr,zi1,zi2,zrs,zis;
  gsl_complex d1,d2,d3,d4,d5,d6;
  
  zr  = fcpt->um[0][0];
  zi1 = fcpt->um[m1][1];
  zi2 = fcpt->um[m2][1];
  zrs = 2*zr;
  zis = zi1+zi2;

  gsl_sf_lngamma_complex_e(1.5-zr,-zi1,&g1,&ag1);
  gsl_sf_lngamma_complex_e(1.5-zr,-zi2,&g2,&ag2);
  gsl_sf_lngamma_complex_e(zrs-1.5,zis,&g3,&ag3);
  gsl_sf_lngamma_complex_e(zr,zi1,&g4,&ag4);
  gsl_sf_lngamma_complex_e(zr,zi2,&g5,&ag5);
  gsl_sf_lngamma_complex_e(3-zrs,-zis,&g6,&ag6);

  d1 = gsl_complex_mul(gsl_complex_polar(exp(g1.val), ag1.val),gsl_complex_polar(exp(g2.val), ag2.val));
  d2 = gsl_complex_mul(d1,gsl_complex_polar(exp(g3.val), ag3.val));
  d3 = gsl_complex_mul(gsl_complex_polar(exp(g4.val), ag4.val),gsl_complex_polar(exp(g5.val), ag5.val));
  d4 = gsl_complex_mul(d3,gsl_complex_polar(exp(g6.val), ag6.val));
  d5 = gsl_complex_div(d2,d4);
  d6 = gsl_complex_mul_real(d5, 0.125*pow(M_PI,-1.5));
  
  return d6;
}

gsl_complex M_A(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_mul_real(n12, 2);
  d1 = gsl_complex_sub (gsl_complex_rect(3,0), d1);
  d2 = gsl_complex_mul_real(n12, 7);
  d2 = gsl_complex_sub (gsl_complex_rect(4,0), d2);
  d1 = gsl_complex_mul(d1,d2);
  
  d2 = gsl_complex_mul (n1, n2);
  d2 = gsl_complex_mul_real(d2, 14);
  
  d1 = gsl_complex_div(d1,d2);
  d1 = gsl_complex_mul(d1,I);
    
  return d1;
}

gsl_complex M_At(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_mul_real(n12, 2);
  d1 = gsl_complex_sub (gsl_complex_rect(3,0), d1);
  d2 = gsl_complex_mul_real(n12, 7);
  d2 = gsl_complex_sub (gsl_complex_rect(8,0), d2);
  d1 = gsl_complex_mul(d1,d2);
  
  d2 = gsl_complex_mul (n1, n2);
  d2 = gsl_complex_mul_real(d2, 14);
  
  d1 = gsl_complex_div(d1,d2);
  d1 = gsl_complex_mul(d1,I);
    
  return d1;
}

gsl_complex M_B(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,d3,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_mul_real(n12, 2);
  d2 = gsl_complex_sub (gsl_complex_rect(3,0), d1);
  d3 = gsl_complex_sub (gsl_complex_rect(1,0), d1);
  d1 = gsl_complex_mul_real(n12, 7);
  d1 = gsl_complex_add (gsl_complex_rect(6,0), d1);
  d1 = gsl_complex_mul(d1, d2);
  d1 = gsl_complex_mul(d1, d3);

  d2 = gsl_complex_add_real(n1, 1);
  d2 = gsl_complex_mul(n1,d2);
  d3 = gsl_complex_add_real(n2, 1);
  d3 = gsl_complex_mul(n2,d3);
  d2 = gsl_complex_mul(d2, d3);
  d2 = gsl_complex_mul_real(d2, -28);

  d1 = gsl_complex_div(d1,d2);
  d1 = gsl_complex_mul(d1,I);
    
  return d1;
}

gsl_complex M_Bt(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,d3,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_mul_real(n12, 2);
  d2 = gsl_complex_sub (gsl_complex_rect(3,0), d1);
  d3 = gsl_complex_sub (gsl_complex_rect(1,0), d1);
  d1 = gsl_complex_mul_real(n12, 7);
  d1 = gsl_complex_add (gsl_complex_rect(-2,0), d1);
  d1 = gsl_complex_mul(d1, d2);
  d1 = gsl_complex_mul(d1, d3);

  d2 = gsl_complex_add_real(n1, 1);
  d2 = gsl_complex_mul(n1,d2);
  d3 = gsl_complex_add_real(n2, 1);
  d3 = gsl_complex_mul(n2,d3);
  d2 = gsl_complex_mul(d2, d3);
  d2 = gsl_complex_mul_real(d2, -28);

  d1 = gsl_complex_div(d1,d2);
  d1 = gsl_complex_mul(d1,I);
    
  return d1;
}

gsl_complex M_C(int m1, int m2)
{
  gsl_complex n1,d1,d2,d3,d4;
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);

  d1 = gsl_complex_add_real(n1,1);
  d2 = gsl_complex_add_real(n1,-1);
  d2 = gsl_complex_mul(n1,d2);
  d3 = gsl_complex_add_real(n1, -2);
  d4 = gsl_complex_add_real(n1, -3);
  d1 = gsl_complex_mul(d1,d2);
  d1 = gsl_complex_mul(d1,d3);
  d1 = gsl_complex_mul(d1,d4);
  d1 = gsl_complex_mul_real(d1, -28.0*M_PI/15.0);
  
  d2 = gsl_complex_mul_real(n1, M_PI);
  d2 = gsl_complex_tan(d2);
  
  d1 = gsl_complex_div(d2,d1);
  
  return d1;
}

gsl_complex M_Ct(int m1, int m2)
{
  gsl_complex n1,d1,d2,d3,d4;
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);

  d1 = gsl_complex_add_real(n1,1);
  d2 = gsl_complex_add_real(n1,-1);
  d2 = gsl_complex_mul(n1,d2);
  d3 = gsl_complex_add_real(n1, -2);
  d4 = gsl_complex_add_real(n1, -3);
  d1 = gsl_complex_mul(d1,d2);
  d1 = gsl_complex_mul(d1,d3);
  d1 = gsl_complex_mul(d1,d4);
  d1 = gsl_complex_mul_real(d1, -28.0*M_PI/9.0);
  
  d2 = gsl_complex_mul_real(n1, M_PI);
  d2 = gsl_complex_tan(d2);
  
  d1 = gsl_complex_div(d2,d1);
  
  return d1;
}

gsl_complex M_D(int m1, int m2)
{
  gsl_complex I;
  
  I = gsl_complex_mul_real(gsl_matrix_complex_get(I_store, m1, m2), 2.0);
  
  return I;
}

gsl_complex M_Ex(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,d3,d4,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_add_real(n1, 1);
  d1 = gsl_complex_mul(n1,d1);
  d2 = gsl_complex_add_real(n2, 1);
  d2 = gsl_complex_mul(n2,d2);
  d1 = gsl_complex_mul(d1, d2);

  d2 = gsl_complex_mul_real(n12, 2);
  d3 = gsl_complex_sub (gsl_complex_rect(3,0), d2);
  d4 = gsl_complex_sub (gsl_complex_rect(1,0), d2);
  d2 = gsl_complex_mul(d3,d4);
  
  d1 = gsl_complex_div(d2,d1);
  d1 = gsl_complex_mul(d1,I);
  
  return d1;
}

gsl_complex M_F(int m1, int m2)
{
  gsl_complex n1,n2,n12,d1,d2,I;
  
  I = gsl_matrix_complex_get(I_store, m1, m2);
  
  n1  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m1][1]);
  n2  = gsl_complex_rect(fcpt->um[0][0], fcpt->um[m2][1]);
  n12 = gsl_complex_rect(2*fcpt->um[0][0], fcpt->um[m1][1]+fcpt->um[m2][1]);

  d1 = gsl_complex_mul(n1, n2);
  
  d2 = gsl_complex_mul_real(n12, 2);
  d2 = gsl_complex_sub (gsl_complex_rect(3,0), d2);
  
  d1 = gsl_complex_div(d2,d1);
  d1 = gsl_complex_mul(d1,I);
  
  return d1;
}

double Pk_fastO(double k, gsl_complex (*func)(int,int))
{
  int i,m1;
  double vec[fcpt->N][2];

  gsl_complex d1,d2,d3,d4,sum;
  sum = gsl_complex_rect(0,0);
  
  for(i=0;i<fcpt->N;i++){
    d1 = gsl_complex_rect(-2*fcpt->um[i][0],-2*fcpt->um[i][1]);
    d2 = gsl_complex_mul_real(d1, log(k));
    d3 = gsl_complex_exp(d2);
    d4 = gsl_complex_mul(gsl_complex_rect(fcpt->cm[i][0],fcpt->cm[i][1]), d3);
    vec[i][0] = GSL_REAL(d4);
    vec[i][1] = GSL_IMAG(d4);
  }

  for(m1=0;m1<fcpt->N;m1++){
    gsl_complex z=(*func)(m1,m1);
    d1 = gsl_complex_mul(gsl_complex_rect(vec[m1][0],vec[m1][1]),z);
    sum = gsl_complex_add(sum,d1);
  }

  return k*k*k*GSL_REAL(sum)*pow((double)fcpt->N,-1)*Pk_lin(k);
}

double Pk_fastT(double k, gsl_complex (*func)(int,int))
{
  int i,m1,m2;
  double vec[fcpt->N][2];

  gsl_complex d1,d2,d3,d4,sum;
  sum = gsl_complex_rect(0,0);
  
  for(i=0;i<fcpt->N;i++){
    d1 = gsl_complex_rect(-2*fcpt->um[i][0],-2*fcpt->um[i][1]);
    d2 = gsl_complex_mul_real(d1, log(k));
    d3 = gsl_complex_exp(d2);
    d4 = gsl_complex_mul(gsl_complex_rect(fcpt->cm[i][0],fcpt->cm[i][1]), d3);
    vec[i][0] = GSL_REAL(d4);
    vec[i][1] = GSL_IMAG(d4);
  }
  
  for(m1=0;m1<fcpt->N;m1++){
    for(m2=0;m2<fcpt->N;m2++){
      gsl_complex z=(*func)(m1,m2);
      d1 = gsl_complex_mul(gsl_complex_rect(vec[m1][0],vec[m1][1]),z);
      d2 = gsl_complex_mul(d1,gsl_complex_rect(vec[m2][0],vec[m2][1]));
      sum = gsl_complex_add(sum,d2);
    }
  }

  return k*k*k*GSL_REAL(sum)*pow((double)fcpt->N,-2);
}

void init_bias(int NFFT)
{
  int i,m1,m2;
  double dlogk;
    
  fcpt = FFTLogPT_init(NFFT,1e-6,1e3,-1.6);
  kkk = (double *)malloc(sizeof(double)*fcpt->N);
  
  I_store = gsl_matrix_complex_alloc(fcpt->N, fcpt->N);
  
  FFTLogPT(fcpt, Pk_lin);
  dlogk = log(fcpt->max/fcpt->min)/(double)fcpt->N;
  
  for(m1=0;m1<fcpt->N;m1++) {
    for(m2=0;m2<fcpt->N;m2++) {
      gsl_matrix_complex_set(I_store, m1 ,m2, Ifunc(m1,m2));
    }
  }

  for(i=0;i<fcpt->N;i++) kkk[i] = fcpt->min*exp((double)i*dlogk);
}

void free_bias(void)
{
    gsl_matrix_complex_free(I_store);
    free(kkk);
}

int comp_bias(int NFFT, double *k, int nk,
	      double *p1,  double *p2, double *p3, double *p4, double *p5, double *p6, double *p7, double *p8, double *p9)
{
  int i;
  double pkt[9][NFFT];
  
  init_bias(NFFT);

  #pragma omp parallel for
  for(i=0;i<fcpt->N;i++) {
    #pragma omp atomic write
    pkt[0][i]=Pk_fastT(kkk[i],&M_A);
    pkt[1][i]=Pk_fastT(kkk[i],&M_At);
    pkt[2][i]=Pk_fastT(kkk[i],&M_B);
    pkt[3][i]=Pk_fastT(kkk[i],&M_Bt);
    pkt[4][i]=Pk_fastO(kkk[i],&M_C);
    pkt[5][i]=Pk_fastO(kkk[i],&M_Ct);
    pkt[6][i]=(Pk_fastT(kkk[i],&M_D)-Pk_fastT(kkk[0],&M_D))*exp(-kkk[i]);
    pkt[7][i]=Pk_fastT(kkk[i],&M_Ex);
    pkt[8][i]=Pk_fastT(kkk[i],&M_F);
    //printf("%le %le %le %le %le %le %le %le %le %le\n",kkk[i],pkA[i],pkAt[i],pkB[i],pkBt[i],pkC[i],pkCt[i],pkD[i],pkE[i],pkF[i]);
  }

  for(i=0;i<9;i++) {
    acc[i] = gsl_interp_accel_alloc();
    spline[i] = gsl_spline_alloc(gsl_interp_cspline,fcpt->N);
    gsl_spline_init(spline[i],kkk,pkt[i],fcpt->N);
  }

  for(i=0;i<nk;i++) {
    p1[i]=gsl_spline_eval(spline[0], k[i], acc[0]);
    p2[i]=gsl_spline_eval(spline[1], k[i], acc[1]);
    p3[i]=gsl_spline_eval(spline[2], k[i], acc[2]);
    p4[i]=gsl_spline_eval(spline[3], k[i], acc[3]);
    p5[i]=gsl_spline_eval(spline[4], k[i], acc[4]);
    p6[i]=gsl_spline_eval(spline[5], k[i], acc[5]);
    p7[i]=gsl_spline_eval(spline[6], k[i], acc[6]);
    p8[i]=gsl_spline_eval(spline[7], k[i], acc[7]);
    p9[i]=gsl_spline_eval(spline[8], k[i], acc[8]);
  }
  
  for(i=0;i<9;i++) {
    gsl_interp_accel_free(acc[i]);
    gsl_spline_free(spline[i]);
  }
  
  free_bias();
    
  return 0;
}
