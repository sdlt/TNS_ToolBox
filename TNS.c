#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <omp.h>
#include "fftwlog.h"

/***! Global variables !***/

NcmFftlog *fcw0;
NcmFftlog *fcw2;
NcmFftlog *fcw4;

static gsl_interp_accel *acc[30];
static gsl_spline *spline[30];
static double kmin, kmax;

/***! P(k) interpolations !***/

double P(double k, int type)
{
  double out;
  
  if (k<kmin) out=gsl_spline_eval(spline[type], kmin, acc[type])*pow(k/kmin,0.9624);
  else if (k>kmax) out=gsl_spline_eval(spline[type], kmax, acc[type])*pow(k/kmax,-2.7);
  else out=gsl_spline_eval(spline[type], k, acc[type]);

  return out;
}

double Pl_2 (double x) {return (3.0*x*x-1.0)*0.5;}
double Pl_4 (double x) {return (35.0*x*x*x*x-30.0*x*x+3.0)*0.125;}

double P_TNS(double mu, void *params)
{
  double Pout=1;
  double mu2,mu4,mu6,mu8;
  double P_gg,P_gt;
  double Leg=1;
  double cA11,cA12,cA22,cA23,cA33,cB111,cB112,cB121,cB122,cB211,cB212,cB221,cB222,cB312,cB321,cB322,cB422;

  double k=((double *)params)[0];
  int l=((double *)params)[1];
  double fg=((double *)params)[2];
  double b1=((double *)params)[3];
  double b2=((double *)params)[4];
  double bg=((double *)params)[5];
  double bt=((double *)params)[6];
  double sigv=((double *)params)[7];
  double alpha_par=((double *)params)[8];
  double alpha_per=((double *)params)[9];
  
  if (l==2) Leg=Pl_2(mu);
  if (l==4) Leg=Pl_4(mu);

  if (fabs(alpha_par*alpha_per-1)>1e-15) {
    double iFe = alpha_per/alpha_par;
    double brac = sqrt(1+mu*mu*(iFe*iFe-1));
    double a0 = pow(alpha_par*alpha_per*alpha_per,1./3.);
    
    Pout/=a0*a0*a0;
    k*=brac/alpha_per;
    mu*=iFe/brac;
  }

  mu2 = mu*mu;
  mu4 = mu2*mu2;
  mu6 = mu4*mu2;
  mu8 = mu4*mu4;

  P_gg = b1*b1*P(k,1)+b1*b2*P(k,2)+2*b1*bg*P(k,4)+2*b1*(bg+2.0/5.0*bt)*P(k,6)+0.25*b2*b2*P(k,8)+bg*bg*P(k,9)+0.5*b2*bg*P(k,10);
  P_gt = b1*P(k,11)+0.5*b2*P(k,3)+bg*P(k,5)+(bg+bt*2.0/5.0)*P(k,7);

  cA11 = pow(b1,3-1)*pow(fg,1) * P(k,13);
  cA12 = pow(b1,3-2)*pow(fg,2) * P(k,14);
  cA22 = pow(b1,3-2)*pow(fg,2) * P(k,15);
  cA23 = pow(b1,3-3)*pow(fg,3) * P(k,16);
  cA33 = pow(b1,3-3)*pow(fg,3) * P(k,17);
  cB111 = pow(b1,4-1-1)*pow(-fg,1+1) * P(k,18);
  cB112 = pow(b1,4-1-2)*pow(-fg,1+2) * P(k,19);
  cB121 = pow(b1,4-2-1)*pow(-fg,2+1) * P(k,20);
  cB122 = pow(b1,4-2-2)*pow(-fg,2+2) * P(k,21);
  cB211 = pow(b1,4-1-1)*pow(-fg,1+1) * P(k,22);
  cB212 = pow(b1,4-1-2)*pow(-fg,1+2) * P(k,23);
  cB221 = pow(b1,4-2-1)*pow(-fg,2+1) * P(k,24);
  cB222 = pow(b1,4-2-2)*pow(-fg,2+2) * P(k,25);
  cB312 = pow(b1,4-1-2)*pow(-fg,1+2) * P(k,26);
  cB321 = pow(b1,4-2-1)*pow(-fg,2+1) * P(k,27);
  cB322 = pow(b1,4-2-2)*pow(-fg,2+2) * P(k,28);
  cB422 = pow(b1,4-2-2)*pow(-fg,2+2) * P(k,29);
  
  Pout *= P_gg
    +mu2*(2*fg*P_gt+cA11+cA12+cB111+cB112+cB121+cB122)
    +mu4*(fg*fg*P(k,12)+cA22+cA23+cB211+cB212+cB221+cB222)
    +mu6*(cA33+cB312+cB321+cB322)
    +mu8*(cB422);
  
  Pout *= 1.0/(1+(k*mu*sigv)*(k*mu*sigv)); // Lorenzian damping
  
  return Pout*Leg;
}

double getPkl(double k, void *params)
{
  double res;
  //double err;
  //size_t ne;
  gsl_function F;
  double par[10];
  gsl_integration_glfixed_table *t=gsl_integration_glfixed_table_alloc(6);
  //gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  par[0]=k;
  par[1]=((double *)params)[0];
  par[2]=((double *)params)[1];
  par[3]=((double *)params)[2];
  par[4]=((double *)params)[3];
  par[5]=((double *)params)[4];
  par[6]=((double *)params)[5];
  par[7]=((double *)params)[6];
  par[8]=((double *)params)[7];
  par[9]=((double *)params)[8];
   
  F.function = &P_TNS;
  F.params = par;
  //gsl_integration_qng(&F, 0, 1, 5e-1, 1e-6, &res, &err, &ne);
  res = gsl_integration_glfixed(&F,0,1,t);
  //gsl_integration_qags(&F, 0, 1, 1e-4, 1e-6, 100, w, &res, &err);

  //gsl_integration_workspace_free(w);
  gsl_integration_glfixed_table_free(t);
  
  return (2.0*par[1]+1.0)*res;
}

double getPklfast(double k, void *params)
{
    int i,nint;
    double res;
    double par[10];
    double dmu;

    par[0]=k;
    par[1]=((double *)params)[0];
    par[2]=((double *)params)[1];
    par[3]=((double *)params)[2];
    par[4]=((double *)params)[3];
    par[5]=((double *)params)[4];
    par[6]=((double *)params)[5];
    par[7]=((double *)params)[6];
    par[8]=((double *)params)[7];
    par[9]=((double *)params)[8];
   
    nint = 400;
    dmu = 1.0/(double)nint;
    
    res = 0;
    #pragma omp parallel for reduction(+:res)
    for (i=0;i<nint;i++) res += P_TNS(0.5*dmu + i*dmu,par);
  
    return (2.0*par[1]+1.0)*res*dmu;
}

double Pl_fft(double k,void *params)
{
  int l=((double *)params)[0];
  int sign=1;
  
  if (l==2 || l==6) sign=-1;
  
  return sign*getPkl(k,params)*pow(k,1.5);
  //return sign*getPklfast(k,params)*pow(k,1.5);
}

int comp_RSD(double *p,double smin,double smax,int ns,double *xi0,double *xi2,double *xi4)
{
  gsl_function F;
  double par[9];

  par[1]=p[0];
  par[2]=p[1];
  par[3]=p[2];
  par[4]=p[3];
  par[5]=p[4];
  par[6]=p[5];
  par[7]=p[6];
  par[8]=p[7];

  F.function=&Pl_fft;

  par[0] = 0;
  F.params = par;
  xi_fftw_vec_lin(&F,xi0,ns,smin,smax,fcw0,0);

  par[0] = 2;
  F.params = par;
  xi_fftw_vec_lin(&F,xi2,ns,smin,smax,fcw2,0);

  par[0] = 4;
  F.params = par;
  xi_fftw_vec_lin(&F,xi4,ns,smin,smax,fcw4,0);

  return 0;
}

int get_prediction(double in_s_min, double in_s_max, int in_n_rsd, double *out, 
		   double in_f, double in_b1, double in_b2, double in_bg, double in_bt, double in_sigv, double in_alpha_par, double in_alpha_per)
{
  int i;
  double xi0[in_n_rsd+1],xi2[in_n_rsd+1],xi4[in_n_rsd+1];
  double par[8];

  par[0]=in_f;
  par[1]=in_b1;
  par[2]=in_b2;
  par[3]=in_bg;
  par[4]=in_bt;
  par[5]=in_sigv;
  par[6]=in_alpha_par;
  par[7]=in_alpha_per;

  comp_RSD(par,in_s_min,in_s_max,in_n_rsd,xi0,xi2,xi4);

  #pragma omp parallel for 
  for (i=1;i<=3*in_n_rsd;i++) {
    if (i<=in_n_rsd) out[i-1]=xi0[i];
    else if (i<=2*in_n_rsd) out[i-1]=xi2[i-in_n_rsd];
    else out[i-1]=xi4[i-2*in_n_rsd];
  }
  
  return 0;
}

int get_prediction_LL(double in_s_min, double in_s_max, int in_n_rsd, double *out, 
		      double in_f, double in_b1, double in_b2, double in_sigv, double in_alpha_par, double in_alpha_per)
{
  int i;
  double xi0[in_n_rsd+1],xi2[in_n_rsd+1],xi4[in_n_rsd+1];
  double par[8];

  par[0]=in_f;
  par[1]=in_b1;
  par[2]=in_b2;
  par[3]=-2.0/7.0*(in_b1-1.0);
  par[4]=11.0/42.0*(in_b1-1.0);
  par[5]=in_sigv;
  par[6]=in_alpha_par;
  par[7]=in_alpha_per;

  comp_RSD(par,in_s_min,in_s_max,in_n_rsd,xi0,xi2,xi4);

  #pragma omp parallel for 
  for (i=1;i<=3*in_n_rsd;i++) {
    if (i<=in_n_rsd) out[i-1]=xi0[i];
    else if (i<=2*in_n_rsd) out[i-1]=xi2[i-in_n_rsd];
    else out[i-1]=xi4[i-2*in_n_rsd];
  }
  
  return 0;
}

double *get_xil_LL(double in_s_min, double in_s_max, int in_n_rsd,
		   double in_f, double in_b1, double in_b2, double in_sigv, double in_alpha_par, double in_alpha_per)
{
  int i;
  double xi0[in_n_rsd+1],xi2[in_n_rsd+1],xi4[in_n_rsd+1];
  double par[8];
  double *out = malloc(sizeof(double)*3*in_n_rsd);

  par[0]=in_f;
  par[1]=in_b1;
  par[2]=in_b2;
  par[3]=-2.0/7.0*(in_b1-1.0);
  par[4]=11.0/42.0*(in_b1-1.0);
  par[5]=in_sigv;
  par[6]=in_alpha_par;
  par[7]=in_alpha_per;

  comp_RSD(par,in_s_min,in_s_max,in_n_rsd,xi0,xi2,xi4);

  #pragma omp parallel for 
  for (i=1;i<=3*in_n_rsd;i++) {
    if (i<=in_n_rsd) out[i-1]=xi0[i];
    else if (i<=2*in_n_rsd) out[i-1]=xi2[i-in_n_rsd];
    else out[i-1]=xi4[i-2*in_n_rsd];
  }
  
  return out;
}

int GetNumLines(char file[])
{
  FILE *fic;
  int n;
  
  fic=fopen(file,"r");
  if (fic) {
    n=0;
    while(!feof(fic)) {
      fscanf(fic,"%*[^\n]\n");
      n++;
    }
    fclose(fic);
    return n;
  } else {
    fprintf(stderr,"Can't read %s file !\n",file);
    return 0;
  }
}

int load_kfuncs(char pkfile[]) 
{
  int i,n,nit;
  FILE *file;
  char filename[BUFSIZ];
  double *k,*pk,*pkn,*pkb1,*pkb2,*pkb3,*pkb4,*pkb5,*pkb6,*pkb7,*pkb8,*pkb9,*pkdt,*pktt,*pkA1,*pkA2,*pkA3,*pkA4,*pkA5,*pkB1,*pkB2,*pkB3,*pkB4,*pkB5,*pkB6,*pkB7,*pkB8,*pkB9,*pkB10,*pkB11,*pkB12;

  //-- Open power spectrum file
  
  strcpy(filename,pkfile);
  strcat(filename,".kfuncs");
  n = GetNumLines(filename);

  //-- Allocate memory

  k = (double *)malloc(sizeof(double)*n);
  pk = (double *)malloc(sizeof(double)*n);
  pkn = (double *)malloc(sizeof(double)*n);
  pkb1 = (double *)malloc(sizeof(double)*n);
  pkb2 = (double *)malloc(sizeof(double)*n);
  pkb3 = (double *)malloc(sizeof(double)*n);
  pkb4 = (double *)malloc(sizeof(double)*n);
  pkb5 = (double *)malloc(sizeof(double)*n);
  pkb6 = (double *)malloc(sizeof(double)*n);
  pkb7 = (double *)malloc(sizeof(double)*n);
  pkb8 = (double *)malloc(sizeof(double)*n);
  pkb9 = (double *)malloc(sizeof(double)*n);
  pkdt = (double *)malloc(sizeof(double)*n);
  pktt = (double *)malloc(sizeof(double)*n);
  pkA1 = (double *)malloc(sizeof(double)*n);
  pkA2 = (double *)malloc(sizeof(double)*n);
  pkA3 = (double *)malloc(sizeof(double)*n);
  pkA4 = (double *)malloc(sizeof(double)*n);
  pkA5 = (double *)malloc(sizeof(double)*n);
  pkB1 = (double *)malloc(sizeof(double)*n);
  pkB2 = (double *)malloc(sizeof(double)*n);
  pkB3 = (double *)malloc(sizeof(double)*n);
  pkB4 = (double *)malloc(sizeof(double)*n);
  pkB5 = (double *)malloc(sizeof(double)*n);
  pkB6 = (double *)malloc(sizeof(double)*n);
  pkB7 = (double *)malloc(sizeof(double)*n);
  pkB8 = (double *)malloc(sizeof(double)*n);
  pkB9 = (double *)malloc(sizeof(double)*n);
  pkB10 = (double *)malloc(sizeof(double)*n);
  pkB11 = (double *)malloc(sizeof(double)*n);
  pkB12 = (double *)malloc(sizeof(double)*n);

  //-- Read data

  file = fopen(filename,"r");
  for (i=0;i<n;i++) fscanf(file,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
			   &k[i], &pk[i], &pkn[i], &pkdt[i], &pktt[i],
			   &pkb1[i], &pkb2[i], &pkb3[i], &pkb4[i], &pkb5[i], &pkb6[i], &pkb7[i], &pkb8[i], &pkb9[i],
			   &pkA1[i], &pkA2[i], &pkA3[i], &pkA4[i], &pkA5[i],
			   &pkB1[i], &pkB2[i], &pkB3[i], &pkB4[i], &pkB5[i], &pkB6[i], &pkB7[i], &pkB8[i], &pkB9[i], &pkB10[i], &pkB11[i], &pkB12[i]);
  fclose(file);

  kmin = k[0];
  kmax = k[n-1];
  printf("> Power spectra read\n");

  //-- Interpolate
  
  nit = 30; 

  for (i=0;i<nit;i++) {
    acc[i] = gsl_interp_accel_alloc();
    spline[i] = gsl_spline_alloc(gsl_interp_cspline,n);
  }
  
  gsl_spline_init(spline[0],k,pk,n);
  gsl_spline_init(spline[1],k,pkn,n);
  gsl_spline_init(spline[2],k,pkb1,n);
  gsl_spline_init(spline[3],k,pkb2,n);
  gsl_spline_init(spline[4],k,pkb3,n);
  gsl_spline_init(spline[5],k,pkb4,n);
  gsl_spline_init(spline[6],k,pkb5,n);
  gsl_spline_init(spline[7],k,pkb6,n);
  gsl_spline_init(spline[8],k,pkb7,n);
  gsl_spline_init(spline[9],k,pkb8,n);
  gsl_spline_init(spline[10],k,pkb9,n);
  gsl_spline_init(spline[11],k,pkdt,n);
  gsl_spline_init(spline[12],k,pktt,n);
  gsl_spline_init(spline[13],k,pkA1,n);
  gsl_spline_init(spline[14],k,pkA2,n);
  gsl_spline_init(spline[15],k,pkA3,n);
  gsl_spline_init(spline[16],k,pkA4,n);
  gsl_spline_init(spline[17],k,pkA5,n);
  gsl_spline_init(spline[18],k,pkB1,n);
  gsl_spline_init(spline[19],k,pkB2,n);
  gsl_spline_init(spline[20],k,pkB3,n);
  gsl_spline_init(spline[21],k,pkB4,n);
  gsl_spline_init(spline[22],k,pkB5,n);
  gsl_spline_init(spline[23],k,pkB6,n);
  gsl_spline_init(spline[24],k,pkB7,n);
  gsl_spline_init(spline[25],k,pkB8,n);
  gsl_spline_init(spline[26],k,pkB9,n);
  gsl_spline_init(spline[27],k,pkB10,n);
  gsl_spline_init(spline[28],k,pkB11,n);
  gsl_spline_init(spline[29],k,pkB12,n);

  free(k);
  free(pk);
  free(pkn);
  free(pkb1);
  free(pkb2);
  free(pkb3);
  free(pkb4);
  free(pkb5);
  free(pkb6);
  free(pkb7);
  free(pkb8);
  free(pkb9);
  free(pkdt);
  free(pktt);
  free(pkA1);
  free(pkA2);
  free(pkA3);
  free(pkA4);
  free(pkA5);
  free(pkB1);
  free(pkB2);
  free(pkB3);
  free(pkB4);
  free(pkB5);
  free(pkB6);
  free(pkB7);
  free(pkB8);
  free(pkB9);
  free(pkB10);
  free(pkB11);
  free(pkB12);

  printf("> Interpolation done\n");
  
  return 0;
}

void free_spline(int nit)
{
  int i;

  for (i=0;i<nit;i++) {
    gsl_interp_accel_free(acc[i]);
    gsl_spline_free(spline[i]);
  }
}

void init_prediction(char pkfile[])
{
  fcw0 = fftwlog_init(512,1e-5,1e3,0.0,0.5);
  fcw2 = fftwlog_init(512,1e-5,1e3,0.0,2.5);
  fcw2 = fftwlog_init(512,1e-5,1e3,0.0,2.5);
  fcw4 = fftwlog_init(512,1e-5,1e3,0.0,4.5);

  load_kfuncs(pkfile);
}

void free_prediction(void)
{
  fftwlog_free_all(fcw0);
  fftwlog_free_all(fcw2);
  fftwlog_free_all(fcw4);

  free_spline(30);
}
  
/***! Main !***/

int main(int argc, char *argv[]) 
{
  int i;
  double xil[120];

  if (argc<2) {
    fprintf(stderr, "Usage: %s [Pk_file.kfuncs]\n", argv[0]);
    exit(1);
  }

  init_prediction(argv[1]);
  //get_prediction(0, 200, 40, xil, 0.5, 1.1, -0.5, 0.1, -0.1, 5, 1.0, 1.0);
  get_prediction_LL(0, 200, 40, xil, 0.5, 1.1, -0.5, 5, 1.0, 1.0);
  free_prediction();
  
  for (i=0;i<40;i++) printf("%le %le %le %le\n",i*5+2.5,xil[i],xil[40+i],xil[80+i]);
  
  return 0;
}
