#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>

double rlog0=-3.0;
double drlog=0.05;
int nr=120;
double r_itp[120];

gsl_interp_accel *accr[3];
gsl_spline *hfr[3];

int init_halofit(double *kk, double *pk, int npk)
{
  int i,index_k,index_r;
  double rmid;
  
  double sum1,sum2,sum3;
  double anorm=0.5*pow(M_PI,-2);

  double arr1[npk],arr2[npk],arr3[npk];
  double arrr1[nr],arrr2[nr],arrr3[nr];
  
  for (i=0;i<nr;i++) r_itp[i]=pow(10,rlog0+i*drlog+drlog*0.5);

  gsl_interp_accel *acc1;
  gsl_spline *hf1;
  gsl_interp_accel *acc2;
  gsl_spline *hf2;
  gsl_interp_accel *acc3;
  gsl_spline *hf3;

  acc1=gsl_interp_accel_alloc();
  acc2=gsl_interp_accel_alloc();
  acc3=gsl_interp_accel_alloc();
  
  hf1=gsl_spline_alloc(gsl_interp_cspline,npk);
  hf2=gsl_spline_alloc(gsl_interp_cspline,npk);
  hf3=gsl_spline_alloc(gsl_interp_cspline,npk);

  accr[0]=gsl_interp_accel_alloc();
  accr[1]=gsl_interp_accel_alloc();
  accr[2]=gsl_interp_accel_alloc();

  hfr[0]=gsl_spline_alloc(gsl_interp_cspline,nr);
  hfr[1]=gsl_spline_alloc(gsl_interp_cspline,nr);
  hfr[2]=gsl_spline_alloc(gsl_interp_cspline,nr);

  for (index_r=0; index_r<nr; index_r++) {
    rmid = r_itp[index_r];
    
    for (index_k=0; index_k<npk; index_k++) {
      double x2 = kk[index_k]*kk[index_k]*rmid*rmid;
      arr1[index_k] = pk[index_k]*pow(kk[index_k],2)*anorm*exp(-x2);
      arr2[index_k] = pk[index_k]*pow(kk[index_k],2)*anorm*2.*x2*exp(-x2);
      arr3[index_k] = pk[index_k]*pow(kk[index_k],2)*anorm*4.*x2*(1.-x2)*exp(-x2);
    }
    gsl_spline_init(hf1,kk,arr1,npk);
    gsl_spline_init(hf2,kk,arr2,npk);
    gsl_spline_init(hf3,kk,arr3,npk);
    
    sum1 = gsl_spline_eval_integ(hf1,kk[0],kk[npk-1],acc1);
    sum2 = gsl_spline_eval_integ(hf2,kk[0],kk[npk-1],acc2);
    sum3 = gsl_spline_eval_integ(hf3,kk[0],kk[npk-1],acc3);
    
    arrr1[index_r] = sqrt(sum1);
    arrr2[index_r] = -sum2/sum1;
    arrr3[index_r] = -sum2*sum2/sum1/sum1 - sum3/sum1;
  }
  
  gsl_spline_init(hfr[0],r_itp,arrr1,nr);
  gsl_spline_init(hfr[1],r_itp,arrr2,nr);
  gsl_spline_init(hfr[2],r_itp,arrr3,nr);

  gsl_interp_accel_free(acc1);
  gsl_interp_accel_free(acc2);
  gsl_interp_accel_free(acc3);

  gsl_spline_free(hf1);
  gsl_spline_free(hf2);
  gsl_spline_free(hf3);
  
  return 0;
}

int comp_halofit(double Omega_m,double Omega_v,double fnu,double Omega0_m,double w0,double *kk,double *pk,double *pkn,int npk)
{
  int index_k,count;
  double pk_lin,pk_quasi,pk_halo,rk;
  double rknl,rneff,rncur;
  double rmid,diff,rmin,rmax;
  
  double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double anorm=0.5*pow(M_PI,-2);

  /** Bisection for knl,neff,ncur **/
  
  rmin=r_itp[0];
  rmax=r_itp[nr-1];
  count = 0;
  do {
    count++;
    rmid = 0.5*(rmax+rmin);
    
    diff=gsl_spline_eval(hfr[0],rmid,accr[0])-1;
    if (diff>0.001) rmin=rmid;
    else if (diff<-0.001) rmax=rmid;

    if(count>100) {
      printf ("Error: number of iterations exceeded\n");
      diff=1e-4;
    }
  } while (fabs(diff) > 0.001);
  
  rknl  = 1./rmid;
  rneff = -3.-gsl_spline_eval(hfr[1],rmid,accr[1]);
  rncur = -gsl_spline_eval(hfr[2],rmid,accr[2]);

  /** Non-linear power spectrum **/
  
  for (index_k=0; index_k<npk; index_k++) {
    
    rk=kk[index_k];
    pk_lin=pk[index_k]*pow(rk,3)*anorm;
    
    gam=0.1971-0.0843*rneff+0.8460*rncur;
    a=1.5222+2.8553*rneff+2.3706*rneff*rneff+0.9903*rneff*rneff*rneff+ 0.2250*rneff*rneff*rneff*rneff-0.6038*rncur+0.1749*Omega_v*(1.+w0);
    a=pow(10,a);
    b=pow(10, (-0.5642+0.5864*rneff+0.5716*rneff*rneff-1.5474*rncur+0.2279*Omega_v*(1.+w0)));
    c=pow(10, 0.3698+2.0404*rneff+0.8161*rneff*rneff+0.5869*rncur);
    xmu=0.;
    xnu=pow(10,5.2105+3.6902*rneff);
    alpha=fabs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
    beta=2.0379-0.7354*rneff+0.3157*pow(rneff,2)+1.2490*pow(rneff,3)+0.3980*pow(rneff,4)-0.1682*rncur + fnu*(1.081 + 0.395*pow(rneff,2));
    
    if(fabs(1-Omega_m)>0.01) { /*then omega evolution */
      f1a=pow(Omega_m,(-0.0732));
      f2a=pow(Omega_m,(-0.1423));
      f3a=pow(Omega_m,(0.0725));
      f1b=pow(Omega_m,(-0.0307));
      f2b=pow(Omega_m,(-0.0585));
      f3b=pow(Omega_m,(0.0743));
      frac=Omega_v/(1.-Omega_m);
      f1=frac*f1b + (1-frac)*f1a;
      f2=frac*f2b + (1-frac)*f2a;
      f3=frac*f3b + (1-frac)*f3a;
    } else {
      f1=1.;
      f2=1.;
      f3=1.;
    }
    
    y=(rk/rknl);
    pk_halo = a*pow(y,f1*3.)/(1.+b*pow(y,f2)+pow(f3*c*y,3.-gam));
    pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*(0.977-18.015*(Omega0_m-0.3)));
    pk_linaa=pk_lin*(1+fnu*47.48*pow(rk,2)/(1+1.5*pow(rk,2)));
    pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);
    
    pkn[index_k]=(pk_halo+pk_quasi)*pow(rk,-3)/anorm;
  }
  
  return 0;
}

void free_halofit(void)
{
  gsl_interp_accel_free(accr[0]);
  gsl_interp_accel_free(accr[1]);
  gsl_interp_accel_free(accr[2]);

  gsl_spline_free(hfr[0]);
  gsl_spline_free(hfr[1]);
  gsl_spline_free(hfr[2]);
}

/*
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

int main(int argc, char *argv[])
{
    int i,n;
    FILE *file;
    double *k,*pk,*pkn;
    double hp[5];
    char filename[BUFSIZ];

    if (argc!=7) {
      fprintf(stderr, "Usage: %s [Pk_file] [Omega_m] [Omega_m] [fnu] [Omega_m0] [w0]\n", argv[0]);
      exit(1);
    }

    sscanf(argv[2], "%lf", &hp[0]);
    sscanf(argv[3], "%lf", &hp[1]);
    sscanf(argv[4], "%lf", &hp[2]);
    sscanf(argv[5], "%lf", &hp[3]);
    sscanf(argv[6], "%lf", &hp[4]);

    n = GetNumLines(argv[1]);
    k = (double *)malloc(sizeof(double)*n);
    pk = (double *)malloc(sizeof(double)*n);
    pkn = (double *)malloc(sizeof(double)*n);
    
    file = fopen((argv[1]),"r");
    for (i=0;i<n;i++) fscanf(file,"%lf %lf",&k[i],&pk[i]);
    fclose(file);
    
    init_halofit(k,pk,n);
    comp_halofit(hp[0],hp[1],hp[2],hp[3],hp[4],k,pk,pkn,n);
    free_halofit();

    strcpy(filename,argv[1]);
    strcat(filename,".halofit");
    file = fopen(filename,"w"); 
    for (i=0;i<n;i++) fprintf(file,"%le %le %le\n",k[i],pk[i],pkn[i]);
    fclose(file);
    
    free(k);
    free(pk);
    free(pkn);
    
    return 0;
}
*/
