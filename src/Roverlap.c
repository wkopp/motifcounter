#include <R.h>
#include "overlap.h"
#include "score2d.h"

#ifdef IN_R
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

#ifdef NOT_USED
void RkeepAlpha(double alpha) {
  Ralpha=alpha;
}

void RkeepDelta(double *delta, double *deltap, int mlen) {
  int k;
  RdeleteDelta();
  #ifdef IN_R
  Rdelta=Calloc(mlen, double);
  Rdeltap=Calloc(mlen, double);
  #else
  Rdelta=calloc(mlen, sizeof(double));
  Rdeltap=calloc(mlen, sizeof(double));
  #endif
  for (k=0;k<mlen;k++) {
    Rdelta[k]=delta[k];
    Rdeltap[k]=deltap[k];
  }
}

void RkeepBeta(double *beta, double *beta3p, double *beta5p, int mlen){
  int k;
  RdeleteBeta();
  #ifdef IN_R
  Rdelta=Calloc(mlen, double);
  Rbeta=Calloc(mlen, double);
  Rbeta3p=Calloc(mlen, double);
  Rbeta5p=Calloc(mlen, double);
  #else
  Rdelta=calloc(mlen, sizeof(double));
  Rbeta=calloc(mlen, sizeof(double));
  Rbeta3p=calloc(mlen, sizeof(double));
  Rbeta5p=calloc(mlen, sizeof(double));
  #endif
  for (k=0;k<mlen;k++) {
    Rbeta[k]=beta[k];
    Rbeta3p[k]=beta3p[k];
    Rbeta5p[k]=beta5p[k];
  }
}
void RdeleteBeta() {
#ifdef IN_R
  if (Rbeta) Free(Rbeta);
  if (Rbeta3p) Free(Rbeta3p);
  if (Rbeta5p) Free(Rbeta5p);
  #else
  if (Rbeta) free(Rbeta);
  if (Rbeta3p) free(Rbeta3p);
  if (Rbeta5p) free(Rbeta5p);
  #endif
  Rbeta=NULL;
  Rbeta3p=NULL;
  Rbeta5p=NULL;
}

void RdeleteDelta() {
#ifdef IN_R
  if (Rdelta) Free(Rdelta);
  if (Rdeltap) Free(Rdeltap);
  #else
  if (Rdelta) free(Rdelta);
  if (Rdeltap) free(Rdeltap);
  #endif
  Rdelta=NULL;
  Rdeltap=NULL;
}
#endif
#endif
