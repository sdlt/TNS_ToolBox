#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <omp.h>

#define KMIN 1e-5
#define KMAX 1e2

/***! Global variables !***/

int NK=896;
double efkmin,efkmax;
double norm = 1;

gsl_interp_accel *acc;
gsl_spline *spla11;
gsl_spline *spla12;
gsl_spline *spla22;
gsl_spline *spla23;
gsl_spline *spla33;
gsl_spline *splb111;
gsl_spline *splb112;
gsl_spline *splb121;
gsl_spline *splb122;
gsl_spline *splb211;
gsl_spline *splb212;
gsl_spline *splb221;
gsl_spline *splb222;
gsl_spline *splb312;
gsl_spline *splb321;
gsl_spline *splb322;
gsl_spline *splb422;

extern double Pk_lin(double k);
extern double Pk_dt(double k);
extern double Pk_tt(double k);

double inPk_lin(double k)
{return Pk_lin(k)/norm;}

double inPk_dt(double k)
{return Pk_dt(k)/norm;}

double inPk_tt(double k)
{return Pk_tt(k)/norm;}

/***! A functions !***/

double A(double r,double x, int ind)
{
  switch(ind) {
  case 11:
    return -r*r*r/7.0*(x+6*x*x*x+r*r*x*(-3.0+10*x*x)+r*(-3.0+x*x-12*x*x*x*x));
  case 12:
    return r*r*r*r/14.0*(x*x-1.0)*(-1.0+7*r*x-6*x*x);
  case 22:
    return r*r*r/14.0*(r*r*x*(13.0-41*x*x)-4*(x+6*x*x*x)+r*(5.0+9*x*x+42*x*x*x*x));
  case 23:
    return r*r*r*r/14.0*(x*x-1.0)*(-1.0+7*r*x-6*x*x);
  case 33:
    return r*r*r/14.0*(1.0-7*r*x+6*x*x)*(-2*x+r*(-1.0+3*x*x));
  default:
    return 0;
  }
}

double At(double r,double x, int ind)
{
  switch(ind) {
  case 11:
    return 1.0/7.0*(x+r-2*r*x*x)*(3*r+7*x-10*r*x*x);
  case 12:
    return r/14.0*(x*x-1.0)*(3*r+7*x-10*r*x*x);
  case 22:
    return 1.0/14.0*(28*x*x+r*x*(25.0-81*x*x)+r*r*(1.0-27*x*x+54*x*x*x*x));
  case 23:
    return r/14.0*(1.0-x*x)*(r-7*x+6*r*x*x);
  case 33:
    return 1.0/14.0*(r-7*x+6*r*x*x)*(-2*x-r+3*r*x*x);
  default:
    return 0;
  }
}

double a(double r, double x, int ind)
{
  switch(ind) {
  case 11:
    return 1.0/7.0*(-7.0*x*x+r*r*r*x*(-3.0+10.0*x*x)+3.0*r*(x+6.0*x*x*x)+r*r*(6.0-19.0*x*x-8.0*x*x*x*x));
  case 12:
    return r/14.0*(-1.0+x*x)*(6.0*r-7.0*(1.0+r*r)*x+8.0*r*x*x);
  case 22:
    return 1.0/14.0*(-28.0*x*x+r*r*r*x*(-13.0+41.0*x*x)+r*x*(11.0+73.0*x*x)-2.0*r*r*(-9.0+31.0*x*x+20.0*x*x*x*x));
  case 23:
    return r/14.0*(-1.0+x*x)*(6.0*r-7.0*(1.0+r*r)*x+8.0*r*x*x);
  case 33:
    return 1.0/14.0*(7.0*x+r*(-6.0+7.0*r*x-8.0*x*x))*(-2.0*x+r*(-1.0+3.0*x*x));
  default:
    return 0;
  }
}

/***! A correction integrals !***/

double corrA_int100 (double x, void * params)
{
  double k=((double *) params)[0];
  double r=((double *) params)[1];
  int ind=((double *) params)[2];
  double l=sqrt(1.0+r*r-2*r*x);
  return A(r,x,ind)*inPk_lin(k)*inPk_lin(k*l)/(l*l*l*l);
}

double corrA_int1 (double r, void * params)
{
  double f;
  double k=((double *) params)[0];
  int ind=((double *) params)[1];
  double prm[3];
  
  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;

  r=exp(r);
  
  if (k<1e-4) return 0;
  prm[0]=k;
  prm[1]=r;
  prm[2]=ind;
  F.function = &corrA_int100;
  F.params = prm;
  gsl_integration_cquad (&F,-1,1,0,1e-3,w,&f,NULL,NULL);
  gsl_integration_cquad_workspace_free(w);
  
  return r*f;
}

double corrA_int200 (double x, void * params)
{
  double k=((double *) params)[0];
  double r=((double *) params)[1];
  int ind=((double *) params)[2];
  double l=sqrt(1.0+r*r-2*r*x);
  return At(r,x,ind)*inPk_lin(k*r)*inPk_lin(k*l)/(l*l*l*l);
}

double corrA_int2 (double r, void * params)
{
  double f;
  double k=((double *) params)[0];
  int ind=((double *) params)[1];
  double prm[3];

  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;

  r=exp(r);
  
  if (k<1e-4) return 0;
  prm[0]=k;
  prm[1]=r;
  prm[2]=ind;
  F.function = &corrA_int200;
  F.params = prm;
  gsl_integration_cquad (&F,-1,1,0,1e-3,w,&f,NULL,NULL);
  gsl_integration_cquad_workspace_free(w);
  
  return r*f;
}

double corrA_int300 (double x, void * params)
{
  double k=((double *) params)[0];
  double r=((double *) params)[1];
  int ind=((double *) params)[2];
  return a(r,x,ind)*inPk_lin(k*r)*inPk_lin(k)/(1.0+r*r-2*r*x);
}

double corrA_int3 (double r, void * params)
{
  double f;
  double k=((double *) params)[0];
  int ind=((double *) params)[1];
  double prm[3];

  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;

  r=exp(r);
  
  if (k<1e-4) return 0;
  prm[0]=k;
  prm[1]=r;
  prm[2]=ind;
  F.function = &corrA_int300;
  F.params = prm;
  gsl_integration_cquad (&F,-1,1,0,1e-3,w,&f,NULL,NULL);
  gsl_integration_cquad_workspace_free(w);
  
  return r*f;
}

double corrA (double k,int m,int n)
{
  double p1,p2,p3,err;
  size_t nevals;
  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;
  double prm[2];

  prm[0]=k;
  prm[1]=m*10+n;

  F.function = &corrA_int1;
  F.params = prm;
  gsl_integration_cquad (&F,log(KMIN),log(KMAX),0.,1e-6,w,&p1,&err,&nevals);

  F.function = &corrA_int2;
  F.params = prm;
  gsl_integration_cquad (&F,log(KMIN),log(KMAX),0.,1e-6,w,&p2,&err,&nevals);
  
  F.function = &corrA_int3;
  F.params = prm;
  gsl_integration_cquad (&F,log(KMIN),log(KMAX),0.,1e-6,w,&p3,&err,&nevals);

  gsl_integration_cquad_workspace_free(w);

  return (p1+p2+p3)/(2*2*M_PI*M_PI)*(k*k*k); // * mu^(2*m)* f^n
}

/***! B functions !***/

double B(double r,double x,int ind)
{
  switch(ind) {
  case 111:
    return r*r/2.0*(x*x-1);
  case 112:
    return 3*r*r/8.0*(x*x-1)*(x*x-1);
  case 121:
    return 3*r*r*r*r/8.0*(x*x-1)*(x*x-1);
  case 122:
    return 5*r*r*r*r/16.0*(x*x-1)*(x*x-1)*(x*x-1);
  case 211:
    return r/2.0*(r+2*x-3*r*x*x);
  case 212:
    return -3.0/4.0*r*(x*x-1)*(-r-2*x+5*r*x*x);
  case 221:
    return 3.0/4.0*r*r*(x*x-1)*(-2+r*r+6*r*x-5*r*r*x*x);
  case 222:
    return -3.0/16.0*r*r*(x*x-1)*(x*x-1)*(6-30*r*x-5*r*r+35*r*r*x*x);
  case 312:
    return r/8.0*(4*x*(3-5*x*x)+r*(3-30*x*x+35*x*x*x*x));
  case 321:
    return r/8.0*(-8*x+r*(-12+36*x*x+12*r*x*(3-5*x*x)+r*r*(3-30*x*x+35*x*x*x*x)));
  case 322:
    return 3.0/16.0*r*(x*x-1)*(-8*x+r*(-12+60*x*x+20*r*x*(3-7*x*x)+5*r*r*(1-14*x*x+21*x*x*x*x)));
  case 422:
    return r/16.0*(8*x*(-3+5*x*x)-6*r*(3-30*x*x+35*x*x*x*x)+6*r*r*x*(15-70*x*x+63*x*x*x*x)+r*r*r*(5-21*x*x*(5-15*x*x+11*x*x*x*x)));
  default:
    return 0;
  }
}

/***! B correction integrals !***/

double corrB_int100 (double x, void * params)
{
  double pk1,pk2;
  double k=((double *) params)[0];
  double r=((double *) params)[1];
  int ind=((double *) params)[2];
  double l=sqrt(1.0+r*r-2*r*x);
  int N=(ind%100)/10;
  int O=(ind%10);

  if (N==1) pk1=inPk_dt(k*l);
  else pk1=inPk_tt(k*l);
  if (O==1) pk2=inPk_dt(k*r);
  else pk2=inPk_tt(k*r);
 
  /*
  pk1=inPk_lin(k*l);
  pk2=inPk_lin(k*r);
  */
  
  return B(r,x,ind)*pk1*pk2*pow(l,-2*N);
}

double corrB_int1 (double r, void * params)
{
  double f;
  double k=((double *) params)[0];
  int ind=((double *) params)[1];
  double prm[3];
  
  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;

  r=exp(r);
  
  if (k<1e-4) return 0;
  prm[0]=k;
  prm[1]=r;
  prm[2]=ind;
  F.function = &corrB_int100;
  F.params = prm;
  gsl_integration_cquad (&F,-1,1,0,1e-3,w,&f,NULL,NULL);
  gsl_integration_cquad_workspace_free(w);
  
  return r*f;
}

double corrB (double k,int m,int n,int o)
{
  double p,err;
  size_t nevals;
  gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (200);
  gsl_function F;
  double prm[2];

  prm[0]=k;
  prm[1]=m*100+n*10+o;

  F.function = &corrB_int1;
  F.params = prm;
  gsl_integration_cquad (&F,log(KMIN),log(KMAX),0.,1e-6,w,&p,&err,&nevals);
  gsl_integration_cquad_workspace_free(w);

  return p/(2*2*M_PI*M_PI)*(k*k*k); // * mu^(2*m) * (-f)^(n+o)
}

/***! Check corrections !***/

void check_corr(void)
{
  FILE *fic;
  int i;
  double dlk = (log10(KMAX)-log10(KMIN))/NK;

  fic=fopen("Acorrection_check.dat","w");
  for (i=0;i<NK;i++) {
    double k=pow(10,log10(KMIN)+i*dlk);
    fprintf(stderr,"%lf\n",k);
    fprintf(fic,"%le %le %le %le %le %le\n",
	    k,
	    corrA(k,1,1),
	    corrA(k,1,2),
	    corrA(k,2,2),
	    corrA(k,2,3),
	    corrA(k,3,3));
  }
  fclose(fic);

  printf("A corrections done.\n");

  fic=fopen("Bcorrection_check.dat","w");
  for (i=0;i<NK;i++) {
    double k=pow(10,log10(KMIN)+i*dlk);
    fprintf(stderr,"%lf\n",k);
    fprintf(fic,"%le %le %le %le %le %le %le %le %le %le %le %le %le\n",
	    k,
	    corrB(k,1,1,1),
	    corrB(k,1,1,2),
	    corrB(k,1,2,1),
	    corrB(k,1,2,2),
	    corrB(k,2,1,1),
	    corrB(k,2,1,2),
	    corrB(k,2,2,1),
	    corrB(k,2,2,2),
	    corrB(k,3,1,2),
	    corrB(k,3,2,1),
	    corrB(k,3,2,2),
	    corrB(k,4,2,2));
  }
  fclose(fic);

  printf("B corrections done.\n");
}

void check_corr_interp(void)
{
  FILE *fic;
  int i;
  double corrAin(double k,int m,int n);
  double corrBin(double k,int m,int n,int o);
  double dlk = (log10(KMAX)-log10(KMIN))/NK;

  fic=fopen("Acorrection_check.dat","w");
  for (i=0;i<NK;i++) {
    double k=pow(10,log10(KMIN)+i*dlk);
    fprintf(fic,"%le %le %le %le %le %le\n",
	    k,
	    corrAin(k,1,1),
	    corrAin(k,1,2),
	    corrAin(k,2,2),
	    corrAin(k,2,3),
	    corrAin(k,3,3));
  }
  fclose(fic);

  printf("A corrections done.\n");

  fic=fopen("Bcorrection_check.dat","w");
  for (i=0;i<NK;i++) {
    double k=pow(10,log10(KMIN)+i*dlk);
    fprintf(fic,"%le %le %le %le %le %le %le %le %le %le %le %le %le\n",
	    k,
	    corrBin(k,1,1,1),
	    corrBin(k,1,1,2),
	    corrBin(k,1,2,1),
	    corrBin(k,1,2,2),
	    corrBin(k,2,1,1),
	    corrBin(k,2,1,2),
	    corrBin(k,2,2,1),
	    corrBin(k,2,2,2),
	    corrBin(k,3,1,2),
	    corrBin(k,3,2,1),
	    corrBin(k,3,2,2),
	    corrBin(k,4,2,2));
  }
  fclose(fic);

  printf("B corrections done.\n");
}

/***! Compute corrections !***/

void comp_TNScorr(double sig8)
{
  int i,nl;
  FILE *fic;
  char file[500];
  double dlk,x[NK];
  double ca11[NK],ca12[NK],ca22[NK],ca23[NK],ca33[NK];
  double cb111[NK],cb112[NK],cb121[NK],cb122[NK],cb211[NK],cb212[NK],cb221[NK],cb222[NK],cb312[NK],cb321[NK],cb322[NK],cb422[NK];
  double lkmin=log10(KMIN);
  double lkmax=log10(KMAX);

  norm = sig8*sig8;
  
  acc = gsl_interp_accel_alloc ();
  
  dlk = (lkmax-lkmin)/NK;
  
  sprintf(file,"Acorrection.dat");
  fic=fopen(file,"r");
  if (fic) {
    nl=0;
    while(!feof(fic)) {
      char ch = fgetc(fic);
      if(ch == '\n') nl++;
    }
    
    if (nl!=NK) printf("Unexpected number of lines (%d/%d) in %s\n",nl,NK,file);
    if (nl<NK) NK=nl;
    else if (nl>NK) exit(1);
    
    rewind(fic);
    
    for (i=0;i<NK;i++) fscanf(fic,"%le %le %le %le %le %le\n",&x[i],&ca11[i],&ca12[i],&ca22[i],&ca23[i],&ca33[i]);
    fclose(fic);
    printf("# %s read.\n",file);
  } else {
    double permax=0;
    
    #pragma omp parallel for
    for (i=0;i<NK;i++) {
      double per,k=pow(10,lkmin+i*dlk);
      
      x[i]=k;
      ca11[i]=corrA(k,1,1)*norm*norm;
      ca12[i]=corrA(k,1,2)*norm*norm;
      ca22[i]=corrA(k,2,2)*norm*norm;
      ca23[i]=corrA(k,2,3)*norm*norm;
      ca33[i]=corrA(k,3,3)*norm*norm;
      
      per=100*(double)(i+1)/NK;
      if (per>permax) {
	#pragma omp atomic read
	permax=per;
      }
      fprintf(stderr,"- A->%.1lf %% done \r",permax);
    }
    
    fic=fopen(file,"w");
    for (i=0;i<NK;i++) fprintf(fic,"%le %le %le %le %le %le\n",x[i],ca11[i],ca12[i],ca22[i],ca23[i],ca33[i]);
    fclose(fic);
    printf("A corrections computed.\n");
  }

  spla11 = gsl_spline_alloc (gsl_interp_cspline, NK);
  spla12 = gsl_spline_alloc (gsl_interp_cspline, NK);
  spla22 = gsl_spline_alloc (gsl_interp_cspline, NK);
  spla23 = gsl_spline_alloc (gsl_interp_cspline, NK);
  spla33 = gsl_spline_alloc (gsl_interp_cspline, NK);
  
  gsl_spline_init (spla11, x, ca11, NK);
  gsl_spline_init (spla12, x, ca12, NK);
  gsl_spline_init (spla22, x, ca22, NK);
  gsl_spline_init (spla23, x, ca23, NK);
  gsl_spline_init (spla33, x, ca33, NK);

  sprintf(file,"Bcorrection.dat");
  fic=fopen(file,"r");
  if (fic) {
    double k;
    for (i=0;i<NK;i++) {
      fscanf(fic,"%le %le %le %le %le %le %le %le %le %le %le %le %le\n",&k,&cb111[i],&cb112[i],&cb121[i],&cb122[i],&cb211[i],&cb212[i],&cb221[i],&cb222[i],&cb312[i],&cb321[i],&cb322[i],&cb422[i]);
      if (fabs(k-x[i])>1e-9) {
	printf("Error in %s\n",file);
	fclose(fic);
	exit(1);
      }
    }
    fclose(fic);
    printf("# %s read.\n",file);
  } else {
    double permax=0;
    
    #pragma omp parallel for
    for (i=0;i<NK;i++) {
      double per,k=x[i];
      
      cb111[i]=corrB(k,1,1,1)*norm*norm;
      cb112[i]=corrB(k,1,1,2)*norm*norm;
      cb121[i]=corrB(k,1,2,1)*norm*norm;
      cb122[i]=corrB(k,1,2,2)*norm*norm;
      cb211[i]=corrB(k,2,1,1)*norm*norm;
      cb212[i]=corrB(k,2,1,2)*norm*norm;
      cb221[i]=corrB(k,2,2,1)*norm*norm;
      cb222[i]=corrB(k,2,2,2)*norm*norm;
      cb312[i]=corrB(k,3,1,2)*norm*norm;
      cb321[i]=corrB(k,3,2,1)*norm*norm;
      cb322[i]=corrB(k,3,2,2)*norm*norm;
      cb422[i]=corrB(k,4,2,2)*norm*norm;

      per=100*(double)(i+1)/NK;
      if (per>permax) {
	#pragma omp atomic read
	permax=per;
      }
      fprintf(stderr,"- B->%.1lf %% done \r",permax);
    }

    fic=fopen(file,"w");
    for (i=0;i<NK;i++) fprintf(fic,"%le %le %le %le %le %le %le %le %le %le %le %le %le\n",x[i],cb111[i],cb112[i],cb121[i],cb122[i],cb211[i],cb212[i],cb221[i],cb222[i],cb312[i],cb321[i],cb322[i],cb422[i]);
    fclose(fic);
    printf("B corrections computed.\n");
  }

  splb111 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb112 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb121 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb122 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb211 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb212 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb221 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb222 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb312 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb321 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb322 = gsl_spline_alloc (gsl_interp_cspline, NK);
  splb422 = gsl_spline_alloc (gsl_interp_cspline, NK);
  
  gsl_spline_init (splb111, x, cb111, NK);
  gsl_spline_init (splb112, x, cb112, NK);
  gsl_spline_init (splb121, x, cb121, NK);
  gsl_spline_init (splb122, x, cb122, NK);
  gsl_spline_init (splb211, x, cb211, NK);
  gsl_spline_init (splb212, x, cb212, NK);
  gsl_spline_init (splb221, x, cb221, NK);
  gsl_spline_init (splb222, x, cb222, NK);
  gsl_spline_init (splb312, x, cb312, NK);
  gsl_spline_init (splb321, x, cb321, NK);
  gsl_spline_init (splb322, x, cb322, NK);
  gsl_spline_init (splb422, x, cb422, NK);

  efkmin = x[0];
  efkmax = x[NK-1];
}

void free_TNScorr(void)
{
  gsl_interp_accel_free (acc);
  gsl_spline_free (spla11);
  gsl_spline_free (spla12);
  gsl_spline_free (spla22);
  gsl_spline_free (spla23);
  gsl_spline_free (spla33);
  gsl_spline_free (splb111);
  gsl_spline_free (splb112);
  gsl_spline_free (splb121);
  gsl_spline_free (splb122);
  gsl_spline_free (splb211);
  gsl_spline_free (splb212);
  gsl_spline_free (splb221);
  gsl_spline_free (splb222);
  gsl_spline_free (splb312);
  gsl_spline_free (splb321);
  gsl_spline_free (splb322);
  gsl_spline_free (splb422);
}

/***! Corrections interpolations !***/

double corrAin(double k,int m,int n)
{
  int ind=m*10+n;

  if (k<efkmin || k>efkmax) return 0;

  switch(ind) {
  case 11:
    return gsl_spline_eval (spla11, k, acc);
  case 12:
    return gsl_spline_eval (spla12, k, acc);
  case 22:
    return gsl_spline_eval (spla22, k, acc);
  case 23:
    return gsl_spline_eval (spla23, k, acc);
  case 33:
    return gsl_spline_eval (spla33, k, acc);
  default:
    return 0;
  }
}

double corrBin(double k,int m,int n,int o)
{
  int ind=m*100+n*10+o;

  if (k<efkmin || k>efkmax) return 0;

  switch(ind) {
  case 111:
    return gsl_spline_eval (splb111, k, acc);
  case 112:
    return gsl_spline_eval (splb112, k, acc);
  case 121:
    return gsl_spline_eval (splb121, k, acc);
  case 122:
    return gsl_spline_eval (splb122, k, acc);
  case 211:
    return gsl_spline_eval (splb211, k, acc);
  case 212:
    return gsl_spline_eval (splb212, k, acc);
  case 221:
    return gsl_spline_eval (splb221, k, acc);
  case 222:
    return gsl_spline_eval (splb222, k, acc);
  case 312:
    return gsl_spline_eval (splb312, k, acc);
  case 321:
    return gsl_spline_eval (splb321, k, acc);
  case 322:
    return gsl_spline_eval (splb322, k, acc);
  case 422:
    return gsl_spline_eval (splb422, k, acc);
  default:
    return 0;
  }
}
