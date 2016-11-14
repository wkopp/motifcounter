#include <R.h>
#include "overlap.h"
#include "score2d.h"

extern DMatrix *Rpwm, *Rcpwm;
extern double *Rtrans, *Rstation, Rgran, Rsiglevel;
extern int Rorder;

void Roverlap(double *alpha, double *beta, double *beta3p, double *beta5p, 
		double *gamma) {

  int i;
  double dx, pvalue;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!beta||!beta3p||!beta5p) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  dx=(double)Rgran;
  pvalue=(double)Rsiglevel;

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

  computeBetas(beta, beta3p,beta5p,gamma,Rpwm->nrow, 0.0);
  *alpha=gamma[0];

}

void RoverlapSingleStranded(double *alpha, double *beta, double *beta3p, double *beta5p, 
		double *gamma) {

  int i;
  double dx, pvalue;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!beta||!beta3p||!beta5p) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }

  dx=(double)Rgran;
  pvalue=(double)Rsiglevel;

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

  computeBetasSingleStranded(beta, gamma,Rpwm->nrow, 0.0);
  *alpha=gamma[0];

}

