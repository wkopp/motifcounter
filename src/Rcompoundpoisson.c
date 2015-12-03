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

void RcompoundpoissonPape_useGamma(double *gamma, 
  double *hitdistribution, int *nseq, int *lseq, int * mhit, int *mclump) {
  int seqlen, i;
  int maxclumpsize, maxhits;
  double lambda;
  double *theta, extention[3];

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!gamma||!hitdistribution||!nseq||!lseq||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

	seqlen=0;
	for (i=0; i<*nseq; i++) {
		seqlen+=lseq[i]-Rpwm->nrow+1;
	}
  //seqlen=(slen[0]-Rpwm->nrow+1)*nos[0];
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
  double *hitdistribution, int *nseq, int *lseq, int * mhit, int *mclump) {
  int seqlen, i;
  int maxclumpsize, maxhits;
  double lambda;
  double *theta, extention[3];
  double *delta, *deltap;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!alpha||!beta||!beta3p||!beta5p||!hitdistribution||!nseq||!lseq||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  //seqlen=(slen[0]-Rpwm->nrow+1)*nos[0];
	seqlen=0;
	for (i=0; i<*nseq; i++) {
		seqlen+=lseq[i]-Rpwm->nrow+1;
	}
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

void Rcompoundpoisson_useBetaSingleStranded(double *alpha, double *beta, 
  double *beta3p, double *beta5p,
  double *hitdistribution, int *nseq, int *lseq, int * mhit, int *mclump) {
  int i, seqlen;
  int maxclumpsize, maxhits;
  double lambda;
  double *theta, extention[3];
  double *delta, *deltap;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!alpha||!beta||!beta3p||!beta5p||!hitdistribution||!nseq||!lseq||!mhit||!mclump) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  //seqlen=(slen[0]-Rpwm->nrow+1)*nos[0];
	seqlen=0;
	for (i=0; i<*nseq; i++) {
		seqlen+=lseq[i]-Rpwm->nrow+1;
	}
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
      alpha[0]/2.0,theta);

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

