#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef IN_R
#include <R.h>
#endif
#include "overlap.h"
#include "forground.h"
#include "background.h"

#define DEBUG
#undef DEBUG
// delta_k=P(X_k'=0,X_k=0, ...., X_1'=0,X_1=0, X'_0|X_0=1)
// deltap_k=P(X_k'=0,X_k=0, ...., X_1'=0,X_1=0 |X'_0=1)
void computeDeltas(double *delta, double *deltap,
 double *beta, double *beta3p, double *beta5p, int mlen) {
  int i,k;
  for (k=0;k<mlen;k++) {
    delta[k]=1.0;
    deltap[k]=1.0;
    for (i=0; i<=k; i++) {
      delta[k]-=(beta[i] + beta3p[i]);
      deltap[k]-=(beta[i] + beta5p[i]);
    }
    #ifdef DEBUG
    #ifdef IN_R
    Rprintf("d%d=%f, d%d'=%f\n",k,delta[k],k,deltap[k]);
    #else
    printf("d%d=%f, d%d'=%f\n",k,delta[k],k,deltap[k]);
    #endif
    #endif
  }
}

// beta_k=P(X_k=1, ...., X_1'=0,X_1=0|X_0=1)
// beta_3pk=P(X_k'=1,X_k=0, ...., X_1'=0,X_1=0|X_0=1)
// beta_5pk=P(X_k=1, ...., X_1'=0,X_1=0|X_0'=1)
void computeBetas(double *beta, double *beta3p, double *beta5p,
  double *gamma, int mlen, double eps) {
  int i,k;
  beta3p[0]=gamma[mlen]-eps;

  //forward-forward
  for (k=1;k<mlen;k++) {
    beta[k]=gamma[k];
    beta3p[k]=gamma[mlen+k];
    beta5p[k]=gamma[mlen*2+k];
    for (i=0; i<k; i++) {
      beta[k]-=(beta[i]*gamma[k-i] +beta3p[i]*gamma[mlen*2+k-i]);
      if (beta[k]<0.0) beta[k]=0;
      beta3p[k]-=(beta[i]*gamma[mlen+k-i] +beta3p[i]*gamma[k-i]);
      if (beta3p[k]<0.0) beta3p[k]=0;
      beta5p[k]-=(beta[i]*gamma[mlen*2+k-i] +beta5p[i]*gamma[k-i]);
      if (beta5p[k]<0.0) beta5p[k]=0;

    }
  }
  #ifdef DEBUG
  #ifdef IN_R
  for (k=0; k<mlen;k++) Rprintf("beta%d=%f\n",k,beta[k]);
  for (k=0; k<mlen;k++) Rprintf("beta3p%d=%f\n",k,beta3p[k]);
  for (k=0; k<mlen;k++) Rprintf("beta5p%d=%f\n",k,beta5p[k]);
  #else
  for (k=0; k<mlen;k++) printf("beta%d=%f\n",k,beta[k]);
  for (k=0; k<mlen;k++) printf("beta3p%d=%f\n",k,beta3p[k]);
  for (k=0; k<mlen;k++) printf("beta5p%d=%f\n",k,beta5p[k]);
  #endif
  #endif
}

void computeBetasSingleStranded(double *beta, double *beta3p, double *beta5p,
  double *gamma, int mlen, double eps) {
  int i,k;
  
  //beta3p[0]=0.0;

  //forward-forward
  for (k=1;k<mlen;k++) {
    beta[k]=gamma[k];
    //beta3p[k]=0.0;
    //beta5p[k]=0.0;
    for (i=0; i<k; i++) {
      beta[k]-=(beta[i]*gamma[k-i]);
      if (beta[k]<0.0) beta[k]=0;
      //beta3p[k]-=(beta[i]*gamma[mlen+k-i] +beta3p[i]*gamma[k-i]);
      //if (beta3p[k]<0.0) beta3p[k]=0;
      //beta5p[k]-=(beta[i]*gamma[mlen*2+k-i] +beta5p[i]*gamma[k-i]);
      //if (beta5p[k]<0.0) beta5p[k]=0;

    }
  }
  #ifdef DEBUG
  #ifdef IN_R
  for (k=0; k<mlen;k++) Rprintf("beta%d=%f\n",k,beta[k]);
  for (k=0; k<mlen;k++) Rprintf("beta3p%d=%f\n",k,beta3p[k]);
  for (k=0; k<mlen;k++) Rprintf("beta5p%d=%f\n",k,beta5p[k]);
  #else
  for (k=0; k<mlen;k++) printf("beta%d=%f\n",k,beta[k]);
  for (k=0; k<mlen;k++) printf("beta3p%d=%f\n",k,beta3p[k]);
  for (k=0; k<mlen;k++) printf("beta5p%d=%f\n",k,beta5p[k]);
  #endif
  #endif
}

