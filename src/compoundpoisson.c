#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef IN_R
#include <R.h>
#endif
#include "score1d.h"
#include "score2d.h"
#include "compoundpoisson.h"
#include "forground.h"
#include "countdist.h"
#include "background.h"
#include "overlap.h"

//#include "countdist.h"
//#include "overlap.h"

#define EPSILON 1e-20
#define DSTRANDED 2
double computePoissonParameter(int seqlen, int mlen, 
   int maxclump, double alpha,double *theta) {
  int i;
  double ec=0.0;
  for (i=0; i<maxclump; i++) {
    ec+=(theta[i*DSTRANDED] + theta[i*DSTRANDED+1])*(double)(i+1);
  }
//  ec=2.15;
  Rprintf("slen=%d, mlen=%d, alpha=%e, ec=%e, lambda=%e\n", seqlen,mlen,alpha,ec, ((double)2*(seqlen-mlen+1)*(alpha))/(ec));
  return ((double)2*(seqlen-mlen+1)*(alpha))/(ec);
}

double *initTheta(int maxclump) {
#ifdef IN_R
  return Calloc(maxclump*DSTRANDED, double);
  #else
  return calloc(maxclump*DSTRANDED, sizeof(double));
  #endif
}

void deleteTheta(double *theta) {
  #ifdef IN_R
  Free(theta);
  #else
  free(theta);
  #endif
}

void computeInitialClump(double *theta, double *gamma, int mlen) {
  int i;
  theta[0]=1-gamma[mlen]+EPSILON;
  theta[1]=1-gamma[mlen]+EPSILON;

  for (i=1; i<mlen; i++) {
    theta[0]*=(1-gamma[i])*(1-gamma[mlen+i]);
    theta[1]*=(1-gamma[i])*(1-gamma[mlen*2+i]);
  }
}

void computeExtentionFactorsPape(double *xi, double *gamma, int mlen) {
  int k, j;
  double xik;
  // xi
  for (k=1; k<mlen; k++) {
    xik=gamma[k]/(1-gamma[k]);
    xik*=(1-gamma[mlen])/(1-gamma[mlen+k]);
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j])/((1-gamma[k+j])*(1-gamma[mlen+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j]);
    }
    #ifdef DEBUG
    #ifdef IN_R
    Rprintf( "xi%d=%f\n", k, xik);
    #else
    printf( "xi%d=%f\n", k, xik);
    #endif
    #endif
    xi[0]+=xik;
  }

  // xi',3'
  for (k=0; k<mlen; k++) {
    if (k==0) {
      xik=(gamma[mlen+k]-EPSILON)/(1-gamma[mlen+k]+EPSILON);
    } else {
      xik=gamma[mlen+k]/(1-gamma[mlen+k]);
    }
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen*2+j])/((1-gamma[k+j])*(1-gamma[mlen+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen*2+j]);
    }
    #ifdef DEBUG
    #ifdef IN_R
    Rprintf( "xi%d'=%f\n", k, xik);
    #else
    printf( "xi%d'=%f\n", k, xik);
    #endif
    #endif
    xi[1]+=xik;
  }

  // xi,5'
  for (k=1; k<mlen; k++) {
    xik=gamma[mlen*2+k]/(1-gamma[mlen*2+k]);
    xik*=(1-gamma[mlen])/(1-gamma[k]);
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j])/((1-gamma[k+j])*(1-gamma[mlen*2+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j]);
    }
    #ifdef DEBUG
    #ifdef IN_R
    Rprintf ("xi%d,5'=%f\n", k,xik);
    #else
    printf ("xi%d,5'=%f\n", k,xik);
    #endif
    #endif
    xi[2]+=xik;
  }
  #ifdef DEBUG
  #ifdef IN_R
  Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #else
  printf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #endif
  #endif
}

void computeCompoundPoissonDistribution(double lambda, 
  int maxhit, int maxclump, double *theta, double *cp) {
  int i, j;
  double p;
  double normalize=0.0;

  #ifdef DEBUG
  #ifdef IN_R
  Rprintf( "lambda=%f\n",lambda);
  #else
  printf( "lambda=%f\n",lambda);
  #endif
  #endif
  cp[0]=exp(-lambda);

  for (i=1;i<maxhit;i++) {
    p=0.0;
    j=i-maxclump+1;
    j=(j>0) ? j : 0;
    for (; j<i;j++) {
      p+=(i-j)*(theta[(i-j-1)*DSTRANDED]+theta[(i-j-1)*DSTRANDED +1])*cp[j];
    }
    cp[i]=lambda/((double)i)*p;
  }
  for (i=0; i<maxhit; i++) normalize+=cp[i];
  for (i=0; i<maxhit; i++) cp[i]/=normalize;
  

}

void computeInitialClumpKopp(double *theta, double *beta3p, 
  double *delta, double *deltap, int mlen) {
  theta[0]=delta[mlen-1]+EPSILON;
  theta[1]=deltap[mlen-1]*(1-beta3p[0]);
}

#define DEBUG
#undef DEBUG
void computeExtentionFactorsKopp(double *xi, 
  double *delta, double *deltap, double *beta, 
  double *beta3p, double *beta5p, int mlen) {
  int k;

  xi[1]=beta3p[0];
  for (k=1;k<mlen;k++) {
    xi[0]+=beta[k];
    xi[1]+=beta3p[k];
    xi[2]+=beta5p[k];
  }
  xi[1]*=deltap[mlen-1]/(delta[mlen-1]+EPSILON);
  xi[2]*=(delta[mlen-1]+EPSILON)/deltap[mlen-1];
  #ifdef DEBUG
  #ifdef IN_R
  Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #else
  printf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #endif
  #endif

}
void computeTheta(int maxclump, double *theta, 
  double *extention, int mlen) {
  int i;
  double total=0.0;

  total=theta[0]+theta[1];
  //Rprintf("%e\n",total);
  for (i=1;i<maxclump; i++) {
    theta[i*DSTRANDED]=
        extention[0]*theta[(i-1)*DSTRANDED] + 
         extention[2]*theta[(i-1)*DSTRANDED+1];
    theta[i*DSTRANDED+1]=
        extention[1]*theta[(i-1)*DSTRANDED] + 
         extention[0]*theta[(i-1)*DSTRANDED+1];

    total+=(theta[i*DSTRANDED]+theta[i*DSTRANDED+1]);
    #ifdef DEBUG
    #ifdef IN_R
  Rprintf("total=%e, theta=%e, theta'=%e\n",total, theta[i*DSTRANDED], 
  theta[i*DSTRANDED+1]);
  #else
  printf("total=%e, theta=%e, theta'=%e\n",total, theta[i*DSTRANDED], 
  theta[i*DSTRANDED+1]);
  #endif
  #endif
  }

  for (i=0;i<maxclump; i++) {
    theta[i*DSTRANDED]/=total;
    theta[i*DSTRANDED+1]/=total;
    //#define DEBUG
    #ifdef DEBUG
    #ifdef IN_R
    Rprintf( "theta%d=%e, theta%d'=%e\n",i+1,
      theta[i*DSTRANDED], i+1,theta[i*DSTRANDED+1]);
    #else
    printf( "theta%d=%e, theta%d'=%e\n",i+1,
      theta[i*DSTRANDED], i+1,theta[i*DSTRANDED+1]);
    #endif
    #endif
  }
}


#ifndef IN_R
#error "obsolete code! use R version or adapt to R version"
void compoundpoisson_correctedpape(char *bgfile, char *pwmfile, 
   char *output_dist, char *slen, char * mhit, 
   char *mclump, char *pv, char *sth, char *gr) {
  DMatrix pwm1, pwm2;
  double *station, *trans;
  int threshold, seqlen;
  int order, i;
  int maxclumpsize, maxhits;
  double dx, pvalue, lambda;
  double *gamma, *theta, extention[3];
  double *hitdistribution;
  FILE *f, *fout=NULL;

  dx=(double)atof(gr);
  if (sth) threshold=(int)roundl(atof(sth)/dx);
  if (pv) pvalue=atof(pv);
  if (slen) seqlen=atoi(slen); else seqlen=1000;
  if (mclump) maxclumpsize=atoi(mclump); else maxclumpsize=15;
  if (mhit) maxhits=atoi(mhit); else maxhits=50;

    f =fopen(bgfile,"rb");
    readBackground(f,&station, &trans, &order);
    fclose(f);

    f =fopen(pwmfile,"rb");
    readMatrix(f,&pwm1);
    getComplementaryMotif(&pwm1,&pwm2);
    fclose(f);

    //fout =fopen(output_quantile,"wb");
    #ifdef IN_R
    gamma=Calloc(pwm1.nrow*4, double);
    #else
    gamma=calloc(pwm1.nrow*4, sizeof(double));
    #endif

    if (pv!=NULL) computeConditionalOverlappingProbabilities(
      &pwm1, &pwm2, station, trans, fout, &pvalue, NULL, &dx, gamma, order);
    else computeConditionalOverlappingProbabilities(
      &pwm1, &pwm2, station, trans, fout, NULL, &threshold, &dx, gamma, order);

    i=0;
    printf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    for (i=1;i<pwm1.nrow; i++) {
      gamma[i]/=gamma[0];
      printf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    }
    for (i=0;i<pwm1.nrow; i++) {
      gamma[pwm1.nrow+i]/=gamma[0];
      printf("p[X%d'=1|X0=1]=%f\n",i,gamma[pwm1.nrow+i]);
    }
    for (i=0;i<pwm1.nrow; i++) {
      gamma[pwm1.nrow*2+i]/=gamma[0];
      printf("p[X%d=1|X0'=1]=%f\n",i,gamma[pwm1.nrow*2+i]);
    }

    memset(extention, 0, 3*sizeof(double));
    computeExtentionFactorsPape(extention, gamma, pwm1.nrow);
    theta=initTheta(maxclumpsize);

    computeInitialClump(theta, gamma,pwm1.nrow);
    computeTheta(maxclumpsize, theta, extention, pwm1.nrow);

    #ifdef IN_R
    hitdistribution=Calloc(maxhits, double);
    #else
    hitdistribution=calloc(maxhits, sizeof(double));
    #endif

    lambda=computePoissonParameter(seqlen, pwm1.nrow, maxclumpsize, gamma[0],theta);
    computeCompoundPoissonDistribution(lambda, maxhits, maxclumpsize, theta, hitdistribution);
    deleteTheta(theta);

    if(output_dist) fout=fopen(output_dist, "w");
    if(output_dist) writeCountDistribution(fout, hitdistribution, maxhits);
    if (output_dist) fclose(fout);

    #ifdef IN_R
    Free(hitdistribution);
    Free(gamma);
    #else
    free(hitdistribution);
    free(gamma);
    #endif
    deleteBackground(station, trans);

    deleteMatrix(&pwm1);
    deleteMatrix(&pwm2);
}

void compoundpoisson_kopp(char *bgfile, char *pwmfile, 
   char *output_dist, char *slen, char * mhit, 
   char *mclump, char *pv, char *sth,char *gr) {

  DMatrix pwm1, pwm2;
  double *station, *trans;
  int threshold, seqlen;
  int order, i;
  int maxclumpsize, maxhits;
  double dx, pvalue, lambda;
  double *gamma, *theta, extention[3], *hitdistribution;
  double *beta, *beta3p, *beta5p, *delta, *deltap;
  FILE *f, *fout=NULL;

  dx=(double)atof(gr);
  if (sth) threshold=(int)roundl(atof(sth)/dx);
  if (pv) pvalue=atof(pv);
  if (slen) seqlen=atoi(slen); else seqlen=1000;
  if (mclump) maxclumpsize=atoi(mclump); else maxclumpsize=15;
  if (mhit) maxhits=atoi(mhit); else maxhits=50;

    f =fopen(bgfile,"rb");
    readBackground(f,&station, &trans, &order);
    fclose(f);

    f =fopen(pwmfile,"rb");
    readMatrix(f,&pwm1);
    getComplementaryMotif(&pwm1,&pwm2);
    fclose(f);

    //fout =fopen(output_quantile,"wb");
    #ifdef IN_R
    gamma=Calloc(pwm1.nrow*4, double);
    #else
    gamma=calloc(pwm1.nrow*4, sizeof(double));
    #endif

    if (pv!=NULL) computeConditionalOverlappingProbabilities(
      &pwm1, &pwm2, station, trans, fout, &pvalue, NULL, &dx, gamma, order);
    else computeConditionalOverlappingProbabilities(
      &pwm1, &pwm2, station, trans, fout, NULL, &threshold, &dx, gamma, order);

    i=0;
      printf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    for (i=1;i<pwm1.nrow; i++) {
      gamma[i]/=gamma[0];
   //   Rprintf"p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    }
    for (i=0;i<pwm1.nrow; i++) {
      gamma[pwm1.nrow+i]/=gamma[0];
   //   Rprintf"p[X%d'=1|X0=1]=%f\n",i,gamma[pwm1.nrow+i]);
    }
    for (i=0;i<pwm1.nrow; i++) {
      gamma[pwm1.nrow*2+i]/=gamma[0];
   //   Rprintf"p[X%d=1|X0'=1]=%f\n",i,gamma[pwm1.nrow*2+i]);
    }
    #ifdef IN_R
    beta=Calloc(pwm1.nrow,double);
    beta3p=Calloc(pwm1.nrow,double);
    beta5p=Calloc(pwm1.nrow,double);
    delta=Calloc(pwm1.nrow,double);
    deltap=Calloc(pwm1.nrow,double);
    #else
    beta=calloc(pwm1.nrow,sizeof(double));
    beta3p=calloc(pwm1.nrow,sizeof(double));
    beta5p=calloc(pwm1.nrow,sizeof(double));
    delta=calloc(pwm1.nrow,sizeof(double));
    deltap=calloc(pwm1.nrow,sizeof(double));
    #endif


    memset(extention, 0, 3*sizeof(double));

    computeBetas(beta, beta3p,beta5p,gamma,pwm1.nrow, 0.0);
    computeDeltas(delta, deltap, beta, beta3p,beta5p,pwm1.nrow);

    computeExtentionFactorsKopp(extention, delta, deltap, 
      beta, beta3p, beta5p, pwm1.nrow);
    //computeExtentionFactorsKopp(extention, gamma, pwm1.nrow);
    theta=initTheta(maxclumpsize);

    computeInitialClumpKopp(theta, gamma,delta, deltap, pwm1.nrow);
    computeTheta(maxclumpsize, theta, extention, pwm1.nrow);
    #ifdef IN_R
    hitdistribution=Calloc(maxhits, double);
    #else
    hitdistribution=calloc(maxhits, sizeof(double));
    #endif

    lambda=computePoissonParameter(seqlen, pwm1.nrow, 
      maxclumpsize, gamma[0],theta);
    computeCompoundPoissonDistribution(lambda, maxhits, 
      maxclumpsize, theta, hitdistribution);
    deleteTheta(theta);

    if(output_dist) fout=fopen(output_dist, "w");
    if(output_dist) writeCountDistribution(fout, hitdistribution, maxhits);
    if (output_dist) fclose(fout);

#ifdef IN_R
    Free(hitdistribution);
    Free(gamma);
    Free(beta);
    Free(beta3p);
    Free(beta5p);
    Free(delta);
    Free(deltap);
    #else
    free(hitdistribution);
    free(gamma);
    free(beta);
    free(beta3p);
    free(beta5p);
    free(delta);
    free(deltap);
    #endif
    deleteBackground(station, trans);

    deleteMatrix(&pwm1);
    deleteMatrix(&pwm2);
}
#endif
