#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <string.h>
#ifdef IN_R
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#endif

#include "background.h"
#include "forground.h"
#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "countdist.h"
#include "posteriorcount.h"
#include "markovchain.h"

extern DMatrix *Rpwm, *Rcpwm;
extern double *Rstation, *Rtrans;

double OverlapHit(int N, double *beta, double *betap) {
  int i;
  double d=1.0, n=1.0;

  if (N<0 || N>=Rpwm->nrow) error("wrong index, i=%d\n", N);

  // beta ... forward hit
  // betap .. reverse hit either 3p or 5p
  // compute denuminator
  for (i=0; i<N; i++) {
    d-=(beta[i]+betap[i]);
  }
  n=(beta[N]);
  //Rprintf("Dn(%d)=%1.4e/%1.4e\n",N,n,d);
  if (d<=0.0) return 0.0;

  return (n/d);
}

double NoOverlapHit(int N, double *beta, double *betap) {
  int i;
  double d=1.0, n=1.0;

  if (N<0 || N>=Rpwm->nrow) error("wrong index, i=%d\n", N);

  // beta ... forward hit
  // betap .. reverse hit either 3p or 5p
  // compute denuminator
  for (i=0; i<N; i++) {
    d-=(beta[i]+betap[i]);
  }
  n=d-(beta[N]+betap[N]);
  //Rprintf("Dn(%d)=%1.4e/%1.4e\n",N,n,d);
  if (d<=0.0) return 0.0;

  return (n/d);
}

#undef DEBUG
#define DEBUG
static double *Rdist=NULL;
//double markovchain(int n, double *par, void *ex) {
  //double *beta, double *beta3p, double *beta5p) {
void markovchain(double *dist, double *a,
 double *beta, double *beta3p, double *beta5p, int slen) {
  int i, k;
  double *post, *prior;
  double alphacond;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }

  // the states are
  // dist[0] ... p(nohit)
  // dist[1] ... p(Hf)
  // dist[2] ... p(Hr)
  // dist[3 ... 3+M-1] ... p(n0), ... , p(nL)
  // dist[3+M, ..., 3+M+M-2] ... p(n1'), ..., p(nL')
  //
  post=Calloc(2*Rpwm->nrow+2, double);
  //prior=calloc(2*Rpwm->nrow+2, sizeof(double));
  prior=dist;
  alphacond=a[0];
  memset(prior, 0, (2*Rpwm->nrow+2)*sizeof(double));
  prior[0]=1-a[0]*2;
  prior[1]=a[0];
  prior[2]=a[0];
  //tmphits[0]=1;

  for (k=0; k<slen; k++) {
    // P(N)
    post[0]=(1-alphacond*(2-beta3p[0]))*(prior[0]+prior[Rpwm->nrow+2] + 
      prior[2*Rpwm->nrow+1]);

    // P(Hf)
    post[1]=alphacond*(prior[0]+prior[Rpwm->nrow+2] +
          prior[2*Rpwm->nrow+1]);

    for (i=1;i<Rpwm->nrow;i++) {
     post[1]+=OverlapHit(i, beta, beta3p)*prior[3+i-1];
    }
    for (i=2;i<Rpwm->nrow;i++) {
     post[1]+=OverlapHit(i, beta5p, beta)*prior[Rpwm->nrow+3+i-2];
    }
    post[1]+=beta5p[1]*prior[2];


    // P(Hr)
    post[2]=alphacond*(1-beta3p[0])*(prior[0]+prior[Rpwm->nrow+2] +
          prior[2*Rpwm->nrow+1]);
    for (i=2;i<Rpwm->nrow;i++) {
     post[2]+=OverlapHit(i, beta, beta5p)*prior[Rpwm->nrow+3+i-2];
    }
    for (i=1;i<Rpwm->nrow;i++) {
     post[2]+=OverlapHit(i, beta3p, beta)*prior[3+i-1];
    }
    // should i switch this line
    post[2]+=beta3p[0]*prior[1];
    //post[2]+=beta3p[0]*post[1];
    post[2]+=beta[1]*prior[2];

    // P(n0)
    post[3]=NoOverlapHit(0, beta, beta3p)*prior[1];
    for (i=1;i<Rpwm->nrow;i++) {
      post[3+i]=NoOverlapHit(i, beta,beta3p)*prior[3+i-1];
    }
    // P(n1')
    post[3+Rpwm->nrow]=NoOverlapHit(1, beta, beta5p)*prior[2];
    for (i=2;i<Rpwm->nrow;i++) {
      post[Rpwm->nrow+3+i-1]=NoOverlapHit(i, beta,beta5p)*prior[Rpwm->nrow+3+i-2];
    }
    //Rprintf("\n");
    memcpy(prior, post, (2*Rpwm->nrow+2)*sizeof(double));
    memset(post, 0, (2*Rpwm->nrow+2)*sizeof(double));
  }

  Free(post);
}

void dmc(int n, double *alphacond, double *gradient, void *ex) {

  double *beta, *beta3p, *beta5p;
  double val;
  double *extra=(double*)ex;
  double epsilon;
  double pa,ma;
  int len;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }

  //alpha=par[0];
  //alpha_expected=extra[0];
  //Rprintf("dmc: alpha_in=%e, alpha_exp=%e\n", alpha, alpha_expected);
  beta=&extra[1];
  beta3p=&extra[Rpwm->nrow+1];
  beta5p=&extra[2*Rpwm->nrow+1];
  len=(int) extra[3*Rpwm->nrow+1];

  if (!Rdist) {
  #ifdef IN_R
    Rdist=Calloc(2*Rpwm->nrow+2, double);
    #else
    Rdist=calloc(2*Rpwm->nrow+2, sizeof(double));
    #endif
  }

//	palpha=exp(log(alphacond[0])+.01);
//	malpha=exp(log(alphacond[0])-.01);
	epsilon=alphacond[0]/1000;
	pa=*alphacond + epsilon;
	ma=*alphacond - epsilon;
  markovchain(Rdist, &pa, beta, beta3p, beta5p, len);

  // determine the gradient numerically
  val=Rdist[1]+Rdist[2];
  markovchain(Rdist, &ma, beta, beta3p, beta5p, len);

  val-=(Rdist[1]+Rdist[2]);
  val/=2*epsilon;

  markovchain(Rdist, alphacond, beta, beta3p, beta5p, len);

  *gradient=-2*(2*extra[0]-Rdist[1]-Rdist[2])*val;
//  Rprintf("grad=%e\n",*gradient);
  //*gradient=val;
}
double minmc(int n, double *alpha, void *ex) {

  double *extra=(double*)ex;
  double *beta, *beta3p, *beta5p;
  int len;
 // alpha=par[0];
 // alpha_expected=(double)ex[0];
  beta=&extra[1];
  beta3p=&extra[Rpwm->nrow+1];
  beta5p=&extra[2*Rpwm->nrow+1];
  len=(int) extra[3*Rpwm->nrow+1];

  if (!Rdist) {
    #ifdef IN_R
    Rdist=Calloc(2*Rpwm->nrow+2, double);
    #else
    Rdist=calloc(2*Rpwm->nrow+2, sizeof(double));
    #endif
  }

  markovchain(Rdist, alpha, beta, beta3p, beta5p,len);

  //Rprintf("minmc: alpha_in=%e, alpha_exp=%e\n", alpha[0], extra[0]);
  //Rprintf("Diff=%e, d^2=%e\n", 2*extra[0]-Rdist[1]-Rdist[2],
//  Rprintf("alpha=%e, alpha'=%e, target=%e\n",2*extra[0], Rdist[1]+Rdist[2],
//  		R_pow_di(extra[0]-Rdist[1], 2));
  return R_pow_di(2*extra[0]-Rdist[1]-Rdist[2], 2);
}

SEXP getMarkovProb(SEXP niter, 
		SEXP alpha, SEXP beta, SEXP beta3p, SEXP beta5p) {
	int i, n;
	int _niter=500;
	double _alpha;
	double *_beta,*_beta3p,*_beta5p;
	double *_rm;
	int nrow;
  SEXP retMatrix;

  //alpha = PROTECT(REAL(alpha));
  _niter=INTEGER(niter)[0];
  _alpha=REAL(alpha)[0];
  beta = PROTECT(coerceVector(beta,REALSXP));
  beta3p = PROTECT(coerceVector(beta3p,REALSXP));
  beta5p = PROTECT(coerceVector(beta5p,REALSXP));

	_beta=REAL(beta);
	_beta3p=REAL(beta3p);
	_beta5p=REAL(beta5p);

	nrow=2*Rpwm->nrow+2;
  PROTECT(retMatrix=allocMatrix(REALSXP,nrow, _niter));

	_rm=REAL(retMatrix);

	if (!Rdist) {
		Rdist=Calloc(nrow, double);
	}
	for (i=0;i<_niter;i++) {
	  markovchain(Rdist, &_alpha, _beta, _beta3p, _beta5p,i);
	  for (n=0;n<nrow;n++) {
	    _rm[i*nrow+n]=Rdist[n];
	  }
	}
  UNPROTECT(4);
  removeDist();
  return retMatrix;
}

void removeDist() {
#ifdef IN_R
  if(Rdist) Free(Rdist);
  #else
  if(Rdist) free(Rdist);
  #endif
  Rdist=NULL;
}

// sampling from the Markov model for debugging purposes:
void sampling_markovchain(double *a,
 double *beta, double *beta3p, double *beta5p, int slen, int nos, int nperm) {
  int i;
  int *cnt, perm,seq,pos;
  int state=0;
  cnt=Calloc(nperm, int);
  double mcnt, vcnt, val, pa;

  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  // state=0 is no hit
  // state=1 is forward strand hit
  // state=2 is reverse strand hit
  GetRNGstate();

  Rprintf("MC-samping:\n");
  Rprintf("alpha=%e, slen=%d, nos=%d, nperm=%d\n", a[0],slen,nos,nperm);

	for (perm=0; perm<nperm; perm++) {
		cnt[perm]=0;
    for (seq=0;seq<nos; seq++) {
      for (pos=0; pos<slen; ) {
        val=unif_rand();
      	if (state==0) {
      		if (val<=1-a[0]*(2-beta3p[0])) {
      			// state remains the same
      		  pos++;
      		} else if (1-a[0]*(2-beta3p[0])< val && val <= 1-a[0]) {
          //Rprintf("%e\n",val);
      			state=2;
      		} else {
          //Rprintf("%e\n",val);
      			state=1;
					}
      	} else if (state==1) {
      		cnt[perm]++;
      		pa=0.0;
      		for (i=0; i<Rpwm->nrow && pa<val; i++) {
      			pa+=beta[i];
      		  if (val<=pa) {
      		  	pos+=i;
      	  		state=1;
      	  		break;
      	  	}
      		}
      		for (i=0; i<Rpwm->nrow && pa<val; i++) {
      			pa+=beta3p[i];
      		  if (val<=pa) {
      		  	pos+=i;
      	  		state=2;
      	  		break;
      		  }
      		}
      		if (pa<val) {
      		  pos+=Rpwm->nrow;
      		  state=0;
      		}
      	} else if (state==2) {
      		cnt[perm]++;
      		pa=0.0;
      		for (i=0; i<Rpwm->nrow && pa<val; i++) {
      			pa+=beta[i];
      		  if (val<=pa) {
      		  	pos+=i;
      	  		state=2;
      	  		break;
      	  	}
      		}
      		for (i=0; i<Rpwm->nrow && pa<val; i++) {
      			pa+=beta5p[i];
      		  if (val<=pa) {
      	  		pos+=i;
      	  		state=1;
      	  		break;
      	  	}
      		}
      		if(pa<val) {
      		  pos+=Rpwm->nrow;
      		  state=0;
      		}
      	}
      }
    }
    //Rprintf("cnt=%d\n",cnt[perm]);
	}
	mcnt=0.0;
	for (i=0; i<nperm; i++) {
		mcnt+=(double)cnt[i];
	}
  //Rprintf("mean(count)=%e\n", mcnt);
	mcnt/=(double)nperm;
	vcnt=0.0;
	for (i=0; i<nperm; i++) {
		vcnt+=(double) ((double)cnt[i]-mcnt)*((double)cnt[i]-mcnt)/((double)(nperm-1));
	}
  // the states are
  // dist[0] ... p(nohit)
  // dist[1] ... p(Hf)
  // dist[2] ... p(Hr)
  // dist[3 ... 3+M-1] ... p(n0), ... , p(nL)
  // dist[3+M, ..., 3+M+M-2] ... p(n1'), ..., p(nL')
  //
  Rprintf("mean(count)=%e\n", mcnt);
  Rprintf("var(count)=%e\n",vcnt);

	PutRNGstate();
}

SEXP sample_mc( 
		SEXP alpha, SEXP beta, SEXP beta3p, SEXP beta5p, SEXP seqlen, SEXP nos, SEXP nperm) {

	double _alpha;
	double *_beta,*_beta3p,*_beta5p;
	int _slen,_nos,_nperm;

  //alpha = PROTECT(REAL(alpha));
  _slen=INTEGER(seqlen)[0];
  _nos=INTEGER(nos)[0];
  _nperm=INTEGER(nperm)[0];
  _alpha=REAL(alpha)[0];
  beta = PROTECT(coerceVector(beta,REALSXP));
  beta3p = PROTECT(coerceVector(beta3p,REALSXP));
  beta5p = PROTECT(coerceVector(beta5p,REALSXP));

	_beta=REAL(beta);
	_beta3p=REAL(beta3p);
	_beta5p=REAL(beta5p);

  sampling_markovchain(&_alpha, _beta, _beta3p, _beta5p, _slen, _nos, _nperm);

  UNPROTECT(3);
	return R_NilValue;
}
