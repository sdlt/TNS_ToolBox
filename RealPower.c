#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

/***! Global variables !***/

static gsl_interp_accel *acc[3];
static gsl_spline *spline[3];
static double kmin, kmax;

extern int comp_bias(int NFFT, double *k, int nk, double *p1,  double *p2, double *p3, double *p4, double *p5, double *p6, double *p7, double *p8, double *p9);
extern int init_halofit(double *kk, double *pk, int npk);
extern int comp_halofit(double Omega_m,double Omega_v,double fnu,double Omega0_m,double w0,double *kk,double *pk,double *pkn,int npk);
extern void free_halofit(void);
extern void comp_TNScorr(double sig8);
extern double corrAin(double k,int m,int n);
extern double corrBin(double k,int m,int n,int o);
extern void free_TNScorr(void);

/***! P(k) interpolations !***/

double Pk_lin(double k)
{
  double out;
  
  if (k<kmin) out=gsl_spline_eval(spline[0], kmin, acc[0])*pow(k/kmin,0.9624);
  else if (k>kmax) out=gsl_spline_eval(spline[0], kmax, acc[0])*pow(k/kmax,-2.7);
  else out=gsl_spline_eval(spline[0], k, acc[0]);

  return out;
}

double Pk_dt(double k)
{
  double out;
  
  if (k<kmin) out=gsl_spline_eval(spline[1], kmin, acc[1])*pow(k/kmin,0.9624);
  else if (k>kmax) out=gsl_spline_eval(spline[1], kmax, acc[1])*pow(k/kmax,-2.7);
  else out=gsl_spline_eval(spline[1], k, acc[1]);

  return out;
}

double Pk_tt(double k)
{
  double out;
  
  if (k<kmin) out=gsl_spline_eval(spline[2], kmin, acc[2])*pow(k/kmin,0.9624);
  else if (k>kmax) out=gsl_spline_eval(spline[2], kmax, acc[2])*pow(k/kmax,-2.7);
  else out=gsl_spline_eval(spline[2], k, acc[2]);

  return out;
}

/*-- Bel's models --*/

double P_dt_Bel(double k, double pk, double pkn, double sig8z)
{
  double rcut = -0.017 + 1.49*sig8z*sig8z;
  double b = 0.091+0.702*sig8z*sig8z;
  double p = sqrt(pk*pkn)*exp(-k*rcut-b*k*k*k*k*k*k);
  
  if (p<0) return 0;
  else return p;
}

double P_tt_Bel(double k, double pk, double sig8z)
{
  double a1 = -0.817 + 3.198*sig8z;
  double a2 = 0.877 - 4.191*sig8z;
  double a3 = -1.199 + 4.629*sig8z;
  double p = pk*exp(-k*(a1+a2*k+a3*k*k));
  
  if (p<0) return 0;
  else return p;
}

/*-- Useful functions --*/

double Sigma2Int(double k, void * params)
{
  double r = *(double *) params;
  double x,y,W2;
  
  k=exp(k);
  x=k*r;
  W2=9*(sin(x)-x*cos(x))*(sin(x)-x*cos(x))/(x*x*x*x*x*x);
  y=k*k*k/(2*M_PI*M_PI)*Pk_lin(k)*W2;
  return y;
}

double GetSigma8(void)
{
  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (100);
  gsl_function F;
  double r=8.0,sig,err;
  
  F.function=&Sigma2Int;
  F.params=&r;
  gsl_integration_cquad (&F,log(kmin),log(kmax),0,1e-9,w,&sig,&err,NULL);
  gsl_integration_cquad_workspace_free(w);
  
  return sqrt(sig);
}

void GetCosmo(double Omega_m0, double redshift, double *Omega_m, double *Dz)
{ 
  double Da,D0,Omega_l,Omega_l0,invE2;

  Omega_l0=1-Omega_m0;
  invE2=1/(Omega_m0*(1+redshift)*(1+redshift)*(1+redshift)+Omega_l0);
  *Omega_m=Omega_m0*(1+redshift)*(1+redshift)*(1+redshift)*invE2;
  Omega_l=Omega_l0*invE2;

  D0=2.5*Omega_m0/(pow(Omega_m0,4./7.)-Omega_l0+(1+Omega_m0*0.5)*(1+Omega_l0/70));
  Da=2.5*(*Omega_m)/(pow((*Omega_m),4./7.)-Omega_l+(1+(*Omega_m)*0.5)*(1+Omega_l/70))/(1+redshift);
  *Dz=Da/D0;
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
    fprintf(stderr,"Can't read %s file!\n",file);
    return 0;
  }
}

int GetNumCols(char file[])
{
  FILE *fic;
  int n;
  
  fic=fopen(file,"r");
  if (fic) {
    n=0;
    while(!feof(fic)) {
      fscanf(fic,"%*s");
      n++;
    }
    fclose(fic);
    return n/GetNumLines(file);
  } else {
    fprintf(stderr,"Can't read %s file!\n",file);
    return 0;
  }
}

/***! Power spectra !***/

void GetPowers(char pkfile[], int mode, double hp[5])
{
  int i,n,nc;
  FILE *file;
  double *k,*pk,*pkn,*pkdt,*pktt;
  double *pkb1,*pkb2,*pkb3,*pkb4,*pkb5,*pkb6,*pkb7,*pkb8,*pkb9;
  char filename[BUFSIZ];
  double sigma8z,fac;

  //-- Check input pk file

  n = GetNumLines(pkfile);
  nc = GetNumCols(pkfile);
  if ((mode==0 && nc!=2) || (mode==1 && nc!=3) || (mode==2 && nc!=5)) {
     fprintf(stderr, "Power spectrum file has a wrong format\n");
     exit(1);
  }

  //-- Allocate memory
    
  k = (double *)malloc(sizeof(double)*n);
  pk = (double *)malloc(sizeof(double)*n);
  pkn = (double *)malloc(sizeof(double)*n);
  pkdt = (double *)malloc(sizeof(double)*n);
  pktt = (double *)malloc(sizeof(double)*n);
  pkb1 = (double *)malloc(sizeof(double)*n);
  pkb2 = (double *)malloc(sizeof(double)*n);
  pkb3 = (double *)malloc(sizeof(double)*n);
  pkb4 = (double *)malloc(sizeof(double)*n);
  pkb5 = (double *)malloc(sizeof(double)*n);
  pkb6 = (double *)malloc(sizeof(double)*n);
  pkb7 = (double *)malloc(sizeof(double)*n);
  pkb8 = (double *)malloc(sizeof(double)*n);
  pkb9 = (double *)malloc(sizeof(double)*n);

   //-- Read linear power
  
  file = fopen(pkfile,"r"); 
  for (i=0;i<n;i++) {
    switch(mode) {
    case 0:
      fscanf(file,"%lf %lf",&k[i],&pk[i]);
      break;
    case 1:
      fscanf(file,"%lf %lf %lf",&k[i],&pk[i],&pkn[i]);
      break;
    case 2:
      fscanf(file,"%lf %lf %lf %lf %lf",&k[i],&pk[i],&pkn[i],&pkdt[i],&pktt[i]);
      break;
    default:
      fprintf(stderr, "Unknown option\n");
      exit(1);
    }
  }
  fclose(file);
  printf("> Input power spectra read\n");
  
  kmin = k[0];
  kmax = k[n-1];

  //-- Compute non linear power

  if (mode==0) {
    init_halofit(k,pk,n);
    comp_halofit(hp[0],hp[1],hp[2],hp[3],hp[4],k,pk,pkn,n);
    free_halofit();
    printf("> Halofit computed\n");
  }

  //-- Interpolate linear power
  
  acc[0] = gsl_interp_accel_alloc();
  spline[0] = gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_spline_init(spline[0],k,pk,n);
  fac = 1;

  //-- Compute Pdt, Ptt and interpolate

  sigma8z = GetSigma8();
  printf("> Using sigma8 = %.4lf\n",sigma8z);

  if (mode==0 || mode==1) {
    for (i=0;i<n;i++) {
      pkdt[i]=P_dt_Bel(k[i],pk[i],pkn[i],sigma8z);
      pktt[i]=P_tt_Bel(k[i],pk[i],sigma8z);
    }
  }

  acc[1] = gsl_interp_accel_alloc();
  spline[1] = gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_spline_init(spline[1],k,pkdt,n);

  acc[2] = gsl_interp_accel_alloc();
  spline[2] = gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_spline_init(spline[2],k,pktt,n);
  
  //-- Compute bias terms
  
  comp_bias(512, k, n, pkb1, pkb2, pkb3, pkb4, pkb5, pkb6, pkb7, pkb8, pkb9);
  printf("> Bias terms computed\n");

  //-- Compute TNS corrections
  
  comp_TNScorr(sigma8z);
  printf("> TNS correction terms computed\n");

  //-- Write real power spectra
  
  strcpy(filename,pkfile);
  strcat(filename,".kfuncs");
  file = fopen(filename,"w"); 
  for (i=0;i<n;i++) fprintf(file,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", k[i], pk[i], pkn[i],
			    pkdt[i], pktt[i], pkb1[i], pkb2[i], pkb3[i], pkb4[i], pkb5[i], pkb6[i], pkb7[i], pkb8[i], pkb9[i],
			    fac*corrAin(k[i],1,1),fac*corrAin(k[i],1,2),fac*corrAin(k[i],2,2),fac*corrAin(k[i],2,3),fac*corrAin(k[i],3,3),
			    fac*corrBin(k[i],1,1,1),fac*corrBin(k[i],1,1,2),fac*corrBin(k[i],1,2,1),fac*corrBin(k[i],1,2,2),fac*corrBin(k[i],2,1,1),fac*corrBin(k[i],2,1,2),
			    fac*corrBin(k[i],2,2,1),fac*corrBin(k[i],2,2,2),fac*corrBin(k[i],3,1,2),fac*corrBin(k[i],3,2,1),fac*corrBin(k[i],3,2,2),fac*corrBin(k[i],4,2,2));
  fclose(file);
  printf("> Power spectra written\n");

  //-- Free
  
  free(k);
  free(pk);
  free(pkn);
  free(pkdt);
  free(pktt);
  free(pkb1);
  free(pkb2);
  free(pkb3);
  free(pkb4);
  free(pkb5);
  free(pkb6);
  free(pkb7);
  free(pkb8);
  free(pkb9);

  for (i=0;i<3;i++) {
    gsl_interp_accel_free(acc[i]);
    gsl_spline_free(spline[i]);
  }
}

/***! Main !***/

int main(int argc, char *argv[]) 
{
  //-- Command-line options
  
  int mode = 0;
  int optind;
  double hp[5];

  if (argc<=3) {
    fprintf(stderr, "Usage: %s [-m] [Pk_file]\n", argv[0]);
    exit(1);
  }
  
  for (optind = 1; optind < argc && argv[optind][0] == '-'; optind++) {
    switch (argv[optind][1]) {
    case 'm':

      if (strcmp(argv[optind+1], "HaloBel") == 0) 
	{
	  // Halofit parameters: Omega_m, Omega_v, fnu, Omega0_m, w0;
	  
	  if (optind+6>argc) {
	    fprintf(stderr, "Missing Halofit parameters\n");
	    exit(1);
	  }
	  sscanf(argv[optind+2], "%lf", &hp[0]);
	  sscanf(argv[optind+3], "%lf", &hp[1]);
	  sscanf(argv[optind+4], "%lf", &hp[2]);
	  sscanf(argv[optind+5], "%lf", &hp[3]);
	  sscanf(argv[optind+6], "%lf", &hp[4]);
	  mode = 0;
	  optind+=6;
	} 
      else if (strcmp(argv[optind+1], "PddBel") == 0)
	{
	  mode = 1;
	  optind++;
	}
      else if (strcmp(argv[optind+1], "AllP") == 0)
	{
	  mode = 2;
	  optind++;
	}
      else
	{
	  fprintf(stderr, "Unknown option\n");
	  exit(1);
	}
      
      break;
    default:
      fprintf(stderr, "Usage: %s [-m] [Pk_file]\n", argv[0]);
      exit(1);
    }   
  }
  argv += optind;
  
  GetPowers(argv[0],mode,hp);
  
  return 0;
}
