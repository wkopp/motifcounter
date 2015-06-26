#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "overlap.h"
//#include "matrix.h"
#include "posteriorcount.h"
#include "score2d.h"
#include "markovchain.h"

extern int Rorder;
extern double *Rstation, *Rtrans;
extern DMatrix *Rcpwm, *Rpwm;
extern double Rgran, Rsiglevel;
//extern double *Rbeta, *Rbeta3p, *Rbeta5p, *Rdelta, *Rdeltap;

#ifdef WK
void dquadratic(int n, double *par, double *gr, void *extra) {
  gr[0]=2*par[0];
}
double quadratic(int n, double *par, void *extra) {
  Rprintf("f=%e\n",par[0]*par[0]);
  return par[0]*par[0];
}
#endif
#define DEBUG
#undef DEBUG
void RPosteriorProbability(double *alpha, double *beta, 
  double *beta3p, double *beta5p,
  double *hitdistribution, int *sseqlen,
  int *smaxhits, int *snos, char **sprior) {
  PosteriorCount prob;
  int seqlen;
  int i, maxhits, k, nos;
  int maxsinglehits;
  double *prior, Zpartition;
  double *singlehitdistribution;
  double *delta, *deltap;
  double a0, aN;
  double psucc,dispersion;
  double abstol=1e-30, intol=1e-30;//, pbio=0.0;
  int trace=0, fail,fncount, type=2, gncount;//, restricted_length;
  double *extra=NULL;
  double zold, znew, eold, enew;
  double sum, res;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!beta||!beta3p||!beta5p||!hitdistribution||
  		!sseqlen||!smaxhits||!snos||!sprior) {
    error("parameters are null");
    return;
  }

  seqlen=sseqlen[0] - Rpwm->nrow+1;
  nos=snos[0];
  maxhits=smaxhits[0];
  maxsinglehits=maxhits;

  delta=Calloc(Rpwm->nrow,double);
  deltap=Calloc(Rpwm->nrow,double);
  if (!delta||!deltap) error("failed to allocate betas");

  computeDeltas(delta, deltap, beta, beta3p,beta5p,Rpwm->nrow);

  // correct bias
  #ifdef IN_R
  extra=Calloc(3*Rpwm->nrow+1, double);
  #else
  extra=calloc(3*Rpwm->nrow+1, sizeof(double));
  #endif
  extra[0]=alpha[0];
  a0=alpha[0];
  for (i=0; i<Rpwm->nrow; i++) {
    extra[1+i]=beta[i];
    extra[Rpwm->nrow+1+i]=beta3p[i];
    extra[2*Rpwm->nrow+1+i]=beta5p[i];
  }

  cgmin(1, &a0, &aN, &res, minmc, dmc, &fail, abstol, intol,
       (void*)extra, type, trace, &fncount, &gncount, 100);
  Rprintf("alpha=%e alpha'=%e\n",a0,aN);

  Free(extra);
  removeDist();
#ifdef SHIFT
#endif

  allocPosteriorProbability(&prob, seqlen, Rpwm->nrow, maxsinglehits);
  initPosteriorProbability(&prob, alpha[0], &beta, &beta3p, &beta5p, 
    &delta, &deltap);

  computePosteriorProbability(&prob);

  prior=Calloc(maxsinglehits+1, double);
#ifdef DEBUG
  if (strcmp(sprior[0],"poisson")==0) {
  	// use poisson from R dist
		for (i=0; i<=maxsinglehits; i++) {
	    prior[i]=dpois((double)i, (double)2*alpha[0]*(seqlen-Rpwm->nrow+1), 0);
		}
  } else if (strcmp(sprior[0],"nbinom")==0) {
  	// use negative binomial
  	dispersion=2;
  	psucc=1/(1+2*alpha[0]*(seqlen-Rpwm->nrow+1)*nos/dispersion);
		for (i=0; i<=maxsinglehits; i++) {
	    prior[i]=dnbinom((double)i, dispersion, psucc, 0);
		}
  } else if (strcmp(sprior[0],"uniform")==0) {
  	// use uniform distribution
    for (i=0; i<=maxsinglehits; i++) {
  	  prior[i]=1.0;
	  }
  } else if (strcmp(sprior[0],"truncunif")==0) {
  	// use truncated uniform distribution
    for (i=0; i<=maxsinglehits; i++) {
  	  prior[i]=1.0;
	  }
  }
#else
  for (i=0; i<=maxsinglehits; i++) {
    prior[i]=1.0;
	}
#endif

  singlehitdistribution=Calloc(maxsinglehits+1, double);

#ifdef DEBUG
    Rprintf("omega=%e, alpha=%e\n", prob.omega, prob.alpha);
#endif


    //singlehitdistribution[0]=addomegas(prob.omega, seqlen, 0)*prior[0];
    singlehitdistribution[0]=prob.probzerohits*prior[0];
    Zpartition=singlehitdistribution[0];
    Rprintf("P(X|%d)=%1.3e\n",0,singlehitdistribution[0]);

    enew=0.0;
    eold=0.0;
    zold=singlehitdistribution[0];
    for (k=1; k<=maxsinglehits; k++) {
      finishPosteriorProbability(&prob, singlehitdistribution, k);


      singlehitdistribution[k]*=prior[k];
      Zpartition+=singlehitdistribution[k];

#ifdef DEBUG
      znew=zold+singlehitdistribution[k];
      enew=zold/znew*eold+(double)k*singlehitdistribution[k]/znew;

			if (strcmp(sprior[0],"truncunif")==0) {
      	if (fabs(enew-(double)2*seqlen*alpha[0]) >
        	 fabs(eold-(double)2*seqlen*alpha[0])) {
//Rprintf("--%d",k);
        	singlehitdistribution[k]=0;
        	break;
      	}
			}
      eold=enew;
      zold=znew;
#endif
    }
    for (k=0; k<=maxsinglehits; k++) {
      singlehitdistribution[k]/=Zpartition;
    }

    multipleShortSequenceProbability(singlehitdistribution, hitdistribution, 
     maxsinglehits, maxhits, nos);
    for (k=0,sum=0.0; k<=maxhits; k++) {
      sum+=hitdistribution[k];
    }

    deletePosteriorProbability(&prob);
    #ifdef IN_R
    Free(singlehitdistribution);
    Free(delta);
    Free(deltap);
    Free(prior);
    #else
    free(singlehitdistribution);
    free(delta);
    free(deltap);
    free(prior);
    #endif

    return;
}
