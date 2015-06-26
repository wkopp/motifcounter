#ifdef IN_R
#include <R.h>
#endif

#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "compoundpoisson.h"

extern double *Rstation, *Rtrans;
extern int Rorder;
extern DMatrix *Rpwm, *Rcpwm;
extern double Rsiglevel, Rgran;

#ifdef WKO
void Rcompoundpoisson_correctedpape(
   double *hitdistribution, int *slen, int *nos, int *mhit, int *mclump) {
  int threshold, seqlen;
  int i;
  int maxclumpsize, maxhits;
  double dx, pvalue, lambda;
  double *gamma, *theta, extention[3];

  if (!Rpwm||!Rcpwm||!Rstation) {
    error("load forground and background properly");
    return;
  }
  if (!hitdistribution||!slen||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  dx=(double)Rgran;
  //threshold=(int)roundl(((double)sth[0])/dx);
  pvalue=(double)Rsiglevel;
  seqlen=slen[0];
  maxclumpsize=(double)mclump[0];
  maxhits=(double)mhit[0];

  //if (!Rbeta) {
  #ifdef IN_R
    gamma=Calloc(Rpwm->nrow*4, double);
  #else
    gamma=calloc(Rpwm->nrow*4, double);
  #endif

    if (pvalue>0.0) 
      computeConditionalOverlappingProbabilities(Rpwm, Rcpwm, 
         Rstation, Rtrans, NULL, 
         &pvalue, NULL, &dx, gamma, Rorder);
    else 
      computeConditionalOverlappingProbabilities(Rpwm, Rcpwm, 
        Rstation, Rtrans, NULL, NULL, &threshold, &dx, gamma, Rorder);

    for (i=1;i<Rpwm->nrow; i++) {
      gamma[i]/=gamma[0];
    }
    for (i=0;i<Rpwm->nrow; i++) {
      gamma[Rpwm->nrow+i]/=gamma[0];
    }
    for (i=0;i<Rpwm->nrow; i++) {
      gamma[Rpwm->nrow*2+i]/=gamma[0];
    }
#ifdef DEBUG
    i=0;
      Rprintf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    for (i=1;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    }
    for (i=0;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d'=1|X0=1]=%f\n",i,gamma[Rpwm->nrow+i]);
    }
    for (i=0;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d=1|X0'=1]=%f\n",i,gamma[Rpwm->nrow*2+i]);
    }
    #endif

    memset(extention, 0, 3*sizeof(double));
    computeExtentionFactorsPape(extention, gamma, Rpwm->nrow);
    theta=initTheta(maxclumpsize);

    computeInitialClump(theta, gamma,Rpwm->nrow);
    computeTheta(maxclumpsize, theta, extention, Rpwm->nrow);

    lambda=computePoissonParameter(seqlen, Rpwm->nrow, 
      maxclumpsize, gamma[0],theta);
    computeCompoundPoissonDistribution(lambda, maxhits, 
      maxclumpsize, theta, hitdistribution);
    deleteTheta(theta);

#ifdef IN_R
    Free(gamma);
    #else
    free(gamma);
    #endif
}
#endif

void RcompoundpoissonPape_useGamma(double *gamma, 
  double *hitdistribution, int *slen, int *nos, int * mhit, int *mclump) {
  int seqlen;
  int maxclumpsize, maxhits;
  double lambda;
  double *theta, extention[3];

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!gamma||!hitdistribution||!slen||!nos||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  seqlen=(slen[0]-Rpwm->nrow+1)*nos[0];
  maxclumpsize=(double)mclump[0];
  maxhits=(double)mhit[0];

  memset(extention, 0, 3*sizeof(double));
  computeExtentionFactorsPape(extention, gamma, Rpwm->nrow);
  Rprintf("ext1=%e, ext2=%e, ext3=%d\n", extention[0],extention[1],
  		extention[2]);
  theta=initTheta(maxclumpsize);

  computeInitialClump(theta, gamma,Rpwm->nrow);
  computeTheta(maxclumpsize, theta, extention, Rpwm->nrow);

  lambda=computePoissonParameter(seqlen, Rpwm->nrow, 
      maxclumpsize, gamma[0],theta);
  computeCompoundPoissonDistribution(lambda, maxhits, 
      maxclumpsize, theta, hitdistribution);
  deleteTheta(theta);
}

void Rcompoundpoisson_useBeta(double *alpha, double *beta, 
  double *beta3p, double *beta5p,
  double *hitdistribution, int *slen, int *nos, int * mhit, int *mclump) {
  int seqlen;
  int maxclumpsize, maxhits;
  double lambda;
  double *theta, extention[3];
  double *delta, *deltap;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!alpha||!beta||!beta3p||!beta5p||!hitdistribution||!slen||!nos||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  seqlen=(slen[0]-Rpwm->nrow+1)*nos[0];
  maxclumpsize=(double)mclump[0];
  maxhits=(double)mhit[0];

#ifdef IN_R
  delta=Calloc(Rpwm->nrow,double);
  deltap=Calloc(Rpwm->nrow,double);
#else
  delta=calloc(Rpwm->nrow,sizeof(double));
  deltap=calloc(Rpwm->nrow,sizeof(double));
#endif

  memset(extention, 0, 3*sizeof(double));

  computeDeltas(delta, deltap, beta, beta3p,beta5p,Rpwm->nrow);


  computeExtentionFactorsKopp(extention, delta, deltap, beta, 
      beta3p, beta5p, Rpwm->nrow);
  theta=initTheta(maxclumpsize);

  computeInitialClumpKopp(theta, beta3p,delta, deltap, Rpwm->nrow);
  computeTheta(maxclumpsize, theta, extention, Rpwm->nrow);

  lambda=computePoissonParameter(seqlen, Rpwm->nrow, maxclumpsize, 
      alpha[0],theta);

  computeCompoundPoissonDistribution(lambda, maxhits, maxclumpsize, 
      theta, hitdistribution);


  deleteTheta(theta);
  #ifdef IN_R
  Free(delta);
  Free(deltap);
  #else
  free(delta);
  free(deltap);
  #endif
}

void Rcompoundpoisson_kopp(
   double *hitdistribution, int *slen, int *nos, int * mhit, int *mclump) {

  int seqlen;
  int i;
  int maxclumpsize, maxhits;
  double dx, pvalue, lambda;
  double *gamma, *theta, extention[3];
  double *beta, *beta3p, *beta5p, *delta, *deltap;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!hitdistribution||!slen||!nos||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  dx=(double)Rgran;
  //threshold=(int)roundl(((double)sth[0])/dx);
  pvalue=(double)Rsiglevel;
  //seqlen=slen[0];
  seqlen=slen[0];
  maxclumpsize=(double)mclump[0];
  maxhits=(double)mhit[0];

    //Rprintf("compute overlapping probabilities\n");
  #ifdef IN_R
    gamma=Calloc(Rpwm->nrow*4, double);
  #else
    gamma=calloc(Rpwm->nrow*4, sizeof(double));
  #endif

    computeConditionalOverlappingProbabilities(Rpwm, Rcpwm, 
        Rstation, Rtrans, NULL, &pvalue, NULL, &dx, gamma, Rorder);

    for (i=1;i<Rpwm->nrow; i++) {
      gamma[i]/=gamma[0];
    }
    for (i=0;i<Rpwm->nrow; i++) {
      gamma[Rpwm->nrow+i]/=gamma[0];
    }
    for (i=0;i<Rpwm->nrow; i++) {
      gamma[Rpwm->nrow*2+i]/=gamma[0];
    }
   #ifdef DEBUG
    i=0;
    Rprintf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    for (i=1;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d=1|X0=1]=%f\n",i,gamma[i]);
    }
    for (i=0;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d'=1|X0=1]=%f\n",i,gamma[Rpwm->nrow+i]);
    }
    for (i=0;i<Rpwm->nrow; i++) {
      Rprintf("p[X%d=1|X0'=1]=%f\n",i,gamma[Rpwm->nrow*2+i]);
    }
    #endif
    #ifdef IN_R
    beta=Calloc(Rpwm->nrow,double);
    beta3p=Calloc(Rpwm->nrow,double);
    beta5p=Calloc(Rpwm->nrow,double);
    delta=Calloc(Rpwm->nrow,double);
    deltap=Calloc(Rpwm->nrow,double);
    #else
    beta=calloc(Rpwm->nrow,sizeof(double));
    beta3p=calloc(Rpwm->nrow,sizeof(double));
    beta5p=calloc(Rpwm->nrow,sizeof(double));
    delta=calloc(Rpwm->nrow,sizeof(double));
    deltap=calloc(Rpwm->nrow,sizeof(double));
    #endif

    memset(extention, 0, 3*sizeof(double));

    computeBetas(beta, beta3p,beta5p,gamma,Rpwm->nrow, 0.0);
    computeDeltas(delta, deltap, beta, beta3p,beta5p,Rpwm->nrow);


    computeExtentionFactorsKopp(extention, delta, deltap, beta, 
      beta3p, beta5p, Rpwm->nrow);
    theta=initTheta(maxclumpsize);

    computeInitialClumpKopp(theta, beta3p,delta, deltap, Rpwm->nrow);
    computeTheta(maxclumpsize, theta, extention, Rpwm->nrow);

    lambda=computePoissonParameter(seqlen, Rpwm->nrow, maxclumpsize, 
      gamma[0],theta);

    computeCompoundPoissonDistribution(lambda, maxhits, maxclumpsize, 
      theta, hitdistribution);
    #ifdef DEBUG
    Rprintf("cpdist=\n");
    for (i=0; i<=maxhits; i++) {
      Rprintf("%1.3e ",hitdistribution[i]);
    }
    Rprintf("\n");
    #endif


    deleteTheta(theta);
    #ifdef IN_R
    Free(gamma);
    Free(beta);
    Free(beta3p);
    Free(beta5p);
    Free(delta);
    Free(deltap);
    #else
    free(gamma);
    free(beta);
    free(beta3p);
    free(beta5p);
    free(delta);
    free(deltap);
    #endif


}
