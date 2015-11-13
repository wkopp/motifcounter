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
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#endif

#include "background.h"
#include "forground.h"
#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "countdist.h"
#include "posteriorcount.h"
#include "markovchain.h"


#undef PARALLEL_WK
#define PARALLEL_WK

int allocPosteriorProbability(PosteriorCount *p, int seqlen, int mlen, int maxhits) {
  int i, j;
  p->seqlen=seqlen;
  p->mlen=mlen;
  p->maxhits=maxhits;
  #ifdef IN_R
  Rprintf("slen=%d, mlen=%d, maxhit=%d\n", p->seqlen, p->mlen, p->maxhits);
  p->value=Calloc(maxhits, double**);
  #else
  p->value=calloc(maxhits, sizeof(double**));
  #endif

  //p->alpha=Calloc(seqlen,double);
  //p->omega=Calloc(seqlen,double);

  if (!p->value) { 
  #ifdef IN_R
    error( "allocPosteriorProbability: error with value calloc\n");
    #else
    fprintf(stderr, "allocPosteriorProbability: error with value calloc\n");
    #endif
    return 1;
  }
  for (i=0; i<maxhits; i++) {
  #ifdef IN_R
    p->value[i]=Calloc(seqlen, double*);
    #else
    p->value[i]=calloc(seqlen, sizeof(double*));
    #endif
    if (!p->value[i]) { 
    #ifdef IN_R
      error("allocPosteriorProbability: error with calloc\n");
      #else
      fprintf(stderr, "allocPosteriorProbability: error with calloc\n");
      #endif
      return 1;
    }
    for (j=0; j<seqlen; j++) {
    #ifdef IN_R
      p->value[i][j]=Calloc(2*mlen, double);
      #else
      p->value[i][j]=calloc(2*mlen, sizeof(double));
      #endif
      if (!p->value[i][j]) {
      #ifdef IN_R
        error("allocPosteriorProbability: error with calloc\n");
        #else
        fprintf(stderr, "allocPosteriorProbability: error with calloc\n");
        #endif
        return 1;
      }
    }
  }
  return 0;
}
void deletePosteriorProbability(PosteriorCount *p) {
#ifdef IN_R
  int i, j;
  for (i=0; i<p->maxhits; i++) {
    for (j=0; j<p->seqlen; j++) Free(p->value[i][j]);
    Free(p->value[i]);
  }
  Free(p->value);
  #else
  int i, j;
  for (i=0; i<p->maxhits; i++) {
    for (j=0; j<p->seqlen; j++) free(p->value[i][j]);
    free(p->value[i]);
  }
  free(p->value);
  #endif
}

double addomegas(double *omega, int start, int end) {
  int i;
	double x=1.0;

  for (i=start; i<=end; i++) {
    x*=omega[i];
  }
  return x;
}

void printVector(double **m, int i1, int i2, int len) {
#ifndef IN_R
  int i;
  for (i=0; i<len;i++) {
    printf("%1.2e\t",m[i1][i2+i]);
  }
  printf("\n");
#endif
}

    #define DEBUG
    #undef DEBUG
void initPosteriorProbability(PosteriorCount *p, double alpha, double **beta,
  double **beta3p, double **beta5p, double **delta, double **deltap) {
  double *extra;
  int i, j;
  int m;
  double abstol=1e-30, intol=1e-30;//, pbio=0.0;
  double res;
  int trace=0, fail,fncount, type=2, gncount;//, restricted_length;
  double a0, aN;
  double *_alpha, *_omega;

  p->beta=*beta;
  p->beta3p=*beta3p;
  p->beta5p=*beta5p;
  p->delta=*delta;
  p->deltap=*deltap;

  _alpha=Calloc(p->seqlen, double);
  _omega=Calloc(p->seqlen, double);

	m=(70>p->seqlen) ? p->seqlen : 70; 

  extra=Calloc(3*p->mlen+2, double);
  extra[0]=alpha;
  a0=alpha;

  for (i=0; i<p->mlen; i++) {
    extra[1+i]=p->beta[i];
    extra[p->mlen+1+i]=p->beta3p[i];
    extra[2*p->mlen+1+i]=p->beta5p[i];
  }

	m=1;
	for (i=0; i<m; i++) {

    extra[3*p->mlen+1]=(double)(500);

		cgmin(1, &a0, &aN, &res, minmc, dmc, &fail, abstol, intol,
				(void*)extra, type, trace, &fncount, &gncount, 100);
		Rprintf("alpha=%e alpha'=%e\n",alpha,aN);

		_alpha[i]=aN;
		a0=aN;
  }

  for (i=m;i<p->seqlen; i++) {
  	_alpha[i]=_alpha[m-1];
  }
  for (i=0; i<p->seqlen; i++) {
    _omega[i]=1-2*_alpha[i]+_alpha[i]*p->beta3p[0];
  }
  p->probzerohits=addomegas(_omega,0,p->seqlen-1);

  
  p->alpha=aN;
  p->omega=1-2*p->alpha+p->alpha*p->beta3p[0];

  Free(extra);
  removeDist();

  #ifdef DEBUG
  #ifdef IN_R
  Rprintf("a=%e, b=%e,b3p=%e, b5p=%e, d=%e, dp=%e, o=%e\n",
  alpha, p->beta[1], p->beta3p[0], p->beta5p[1], p->delta[0], p->deltap[0],
  p->omega);
  #endif
  #endif

  // init k=1
  for (j=0; j<p->seqlen; j++) {
    for (i=0; i<=j; i++) {
      if (p->mlen -1<= j-i) {

        p->value[0][j][0]+=
            _alpha[i] * p->delta[p->mlen-1] 
            * addomegas(_omega, 0, i-1)
            * addomegas(_omega, i+p->mlen, j);

        p->value[0][j][p->mlen]+= 
            _alpha[i] * p->delta[0] * p->deltap[p->mlen-1]
            * addomegas(_omega, 0, i-1)
            * addomegas(_omega, i+p->mlen, j);

      } else {
        p->value[0][j][p->mlen-j+i-1]+= _alpha[i] * addomegas(_omega,0, i-1);
        p->value[0][j][2*p->mlen-j+i-1]+= _alpha[i] * 
        	p->delta[0] * addomegas(_omega,0, i-1);
      }
    }

  //#define DEBUG
    #ifdef DEBUG
    #ifndef IN_R
    for (i=0; i<2*p->mlen; i++) {
      printf( "%1.2e\t",p->value[0][j][i]);
    }
    printf( "\n");
    #else
    for (i=0; i<2*p->mlen; i++) {
      Rprintf( "%1.2e\t",p->value[0][j][i]);
    }
    Rprintf( "\n");
    #endif
    #endif
  }
  Free(_alpha);
  Free(_omega);
}

double ffTransProb(PosteriorCount *prob, int n) {
  if (n<0 ) return 0.0;
  else if (n<prob->mlen) return (prob->beta[n]);
  else return prob->alpha;
}

double rrTransProb(PosteriorCount *prob, int n) {
  if (n<0 ) return 0.0;
  else if (n<prob->mlen) return (prob->beta[n]);
  else return prob->alpha*prob->delta[0];
}

double frTransProb(PosteriorCount *prob, int n) {
  if (n<0 ) return 0.0;
  else if (n<prob->mlen) return (prob->beta3p[n]);
  else return prob->alpha*prob->delta[0];
}

double rfTransProb(PosteriorCount *prob, int n) {
  if (n<0 ) return 0.0;
  else if (n<prob->mlen) return (prob->beta5p[n]);
  else return prob->alpha ;
}

double rNonHitStretch(PosteriorCount *prob, int n) {
  if (n<prob->mlen) return 1.0;
  else {
    return prob->deltap[prob->mlen-1] * R_pow_di(prob->omega, n-prob->mlen);
  }
}

double fNonHitStretch(PosteriorCount *prob, int n) {
  if (n<prob->mlen) return 1.0;
  else return prob->delta[prob->mlen-1] * R_pow_di(prob->omega, n-prob->mlen);
}

    #define DEBUG
    #undef DEBUG
void computePosteriorProbability(PosteriorCount *prob) {
  int i, j, k;
  int mc, mp;
  #ifdef DEBUG
  int p;
  #endif
  for (k=1; k<prob->maxhits; k++) {
  #ifndef IN_R
  #ifdef DEBUG
    fprintf(stdout, "k=%d ... \n", k);
  #endif
  #endif
  	//Rprintf("R_pow_di(0.1,0)=%e\n",R_pow_di(0.1,0));
  #ifdef PARALLEL_WK
  #ifdef _OPENMP
   #pragma omp parallel for default(none) shared(k,prob) private(i,j, mc,mp)
  #endif
  #endif
    for (i=0; i<prob->seqlen; i++) {
      R_CheckUserInterrupt();
      for (j=0; j<i; j++) {

        if (prob->mlen <= i-j) mc=0;
        else mc=prob->mlen-i+j;

        for (mp=0; mp<prob->mlen; mp++) {
          // forward hit
          prob->value[k][i][mc]+=
              (prob->value[k-1][j][mp] * ffTransProb(prob, prob->mlen-mp) +
               prob->value[k-1][j][prob->mlen+mp] * rfTransProb(prob, prob->mlen-mp))
                 *fNonHitStretch(prob, i-j);

          // reverse hit
          prob->value[k][i][prob->mlen+mc]+=
                (prob->value[k-1][j][prob->mlen+mp] * rrTransProb(prob, prob->mlen-mp) +
                            prob->value[k-1][j][mp] * frTransProb(prob, prob->mlen-mp))
                            *rNonHitStretch(prob, i-j);


////////////////////////
//
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
          fprintf(stdout, "p(X_[0,%d],a=%d|k=%d)+=(p([0,%d],%d|%d)*%1.1e + p([0,%d]',%d|%d)*%1.1e)*%1.1e=\n ",
           i,mc,k,
           j,mp,k-1,
          ffTransProb(prob, prob->mlen-mp), 
           j,mp,k-1,
          rfTransProb(prob, prob->mlen-mp), fNonHitStretch(prob, i-j));

          fprintf(stdout, "\t (%1.1e*%1.1e+%1.1e*%1.1e)*%1.1e\n",
          prob->value[k-1][j][mp], 
          ffTransProb(prob, prob->mlen-mp), 
          prob->value[k-1][j][prob->mlen+mp], 
          rfTransProb(prob, prob->mlen-mp), fNonHitStretch(prob, i-j));

          fprintf(stdout, "p(X_[0,%d]',a=%d|k=%d)+=(p([0,%d]',%d|%d)*%1.1e + p([0,%d],%d|%d)*%1.1e)*%1.1e=\n ",
           i,mc,k,
           j,mp,k-1,
          rrTransProb(prob, prob->mlen-mp), 
           j,mp,k-1,
          frTransProb(prob, prob->mlen-mp), fNonHitStretch(prob, i-j));

          fprintf(stdout, "\t (%1.1e*%1.1e+%1.1e*%1.1e)*%1.1e\n",
          prob->value[k-1][j][prob->mlen+mp], 
          rrTransProb(prob, prob->mlen-mp), 
          prob->value[k-1][j][mp], 
          frTransProb(prob, prob->mlen-mp), rNonHitStretch(prob, i-j));
#endif
#endif

        }
      }

      for (j=0; j<=i; j++) {
      	R_CheckUserInterrupt();
        if (prob->mlen-1<= i-j) mc=0;
        else mc=prob->mlen-i+j-1;

        prob->value[k][i][prob->mlen+mc]+=prob->value[k-1][j][prob->mlen-1] * 
          frTransProb(prob, 0)*rNonHitStretch(prob, i-j+1);

#define DEBUG
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
        fprintf(stdout, "p(X_[0,%d]',a=%d|k=%d)+=p([0,%d],%d|%d) * %1.1e * %1.1e=\n ",
           i,mc,k,
           j,prob->mlen-1,k-1,
           //prob->value[k-1][j][mp],
          frTransProb(prob, 0), rNonHitStretch(prob, i-j+1));

        fprintf(stdout, "\t %1.1e * %1.1e * %1.1e\n",
          prob->value[k-1][j][prob->mlen-1], 
          frTransProb(prob, 0), rNonHitStretch(prob, i-j+1));
#endif
#endif
      }

#define DEBUG
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
      for(p=0;p<prob->mlen*2;p++) {
        fprintf(stdout, "%1.2e\t",prob->value[k][i][p]);
      }
      fprintf(stdout, "\n");
#endif
#endif
    }

  }
  //printResult(prob);
}

void printResult(PosteriorCount *prob) {
	int k, p;
	for (k=0; k<prob->maxhits; k++) {
    for(p=0;p<prob->mlen*2;p++) {
      Rprintf("%1.2e\t",prob->value[k][prob->seqlen-1][p]);
    }
    Rprintf("\n");
  }
}

void finishPosteriorProbability(PosteriorCount *prob, double *final, int nhits) {
  int m;
  final[nhits]+=prob->value[nhits-1][prob->seqlen-1][0];
  for (m=1; m<prob->mlen; m++) {
    final[nhits]+=prob->value[nhits-1][prob->seqlen-1][m]*(prob->delta[prob->mlen-m-1]);
  }
  final[nhits]+=prob->value[nhits-1][prob->seqlen-1][prob->mlen];
  for (m=1; m<prob->mlen; m++) {
    final[nhits]+=prob->value[nhits-1][prob->seqlen-1][prob->mlen+m]*(prob->deltap[prob->mlen-m-1]);
  }
#ifdef WK
	for (m=0; m<prob->mlen; m++) {
		final[nhits]+=prob->value[nhits-1][prob->seqlen-1][m]*(prob->delta[prob->mlen-m]);
		final[nhits]+=prob->value[nhits-1][prob->seqlen-1][prob->mlen+m]*(prob->deltap[prob->mlen-m]);
	}
#endif
}

#undef DEBUG
void convolute(double *result, double *p1, double *p2, int len) {
  int i, j;
  for(i=0;i<=len;i++) {
    for (j=0;j<=len && i+j<=len; j++) {
      result[i+j]+=p1[i]*p2[j];
    }
  }
}

#ifdef WK
double onehit(double alpha, double *b,
		double *b3, double *b5, double *d, double *dp, double omega,
		int mlen, int slen) {
	int i;
	double nh;
	double prob=0.0;

	for (i=0; i<slen; i++) {
    if (i+mlen>slen) {
      nh= d[slen-i-2];
    } else {
    	nh=d[mlen-1]*R_pow_di(omega,slen-i-mlen);
    }
    prob+=R_pow_di(omega,i)*alpha*nh;
	}
	for (i=0; i<slen; i++) {
    if (i+mlen>slen) {
      nh= dp[slen-i-2];
    } else {
    	nh=dp[mlen-1]*R_pow_di(omega,slen-i-mlen);
    }
    prob+=R_pow_di(omega,i)*alpha*(1-b3[0])*nh;
	}
	return prob;
}

double twohit(double alpha, double *b,
		double *b3, double *b5, double *d, double *dp, double omega,
		int mlen, int slen) {
	int i, j;
	double prob=0.0;

	// forward
	for (i=0; i<slen; i++) {
		//forward
		for (j=i+1; j<slen; j++) {
		}
		// reverse
		for (j=i+1; j<slen; j++) {
		}
	} 
	// reverse
	for (i=0; i<slen; i++) {
		//forward
		for (j=i+1; j<slen; j++) {
		}
		// reverse
		for (j=i+1; j<slen; j++) {
		}
	} 
}

double threehit(double alpha, double *b,
		double *b3, double *b5, double *d, double *dp, double omega,
		int mlen, int slen) {
	int i,j,k;
	// forward
	for (i=0; i<slen; i++) {
		//forward
		for (j=i+1; j<slen; j++) {
			for (k=j+1; k<slen; k++) {
			}
			for (k=j+1; k<slen; k++) {
			}
		}
		// reverse
		for (j=i+1; j<slen; j++) {
			for (k=j+1; k<slen; k++) {
			}
			for (k=j+1; k<slen; k++) {
			}
		}
	} 
	// reverse
	for (i=0; i<slen; i++) {
		//forward
		for (j=i+1; j<slen; j++) {
			for (k=j+1; k<slen; k++) {
			}
			for (k=j+1; k<slen; k++) {
			}
		}
		// reverse
		for (j=i+1; j<slen; j++) {
			for (k=j+1; k<slen; k++) {
			}
			for (k=j+1; k<slen; k++) {
			}
		}
	} 
}
#endif

void computeResultRecursive(double ** part, int nos, int klen) {
  int l1, l2;
#ifdef DEBUG
	int i;
  double sum;
#endif
  l1=nos/2;
  l2=nos-l1;

#ifdef DEBUG
  Rprintf("nos=%d from l1=%d l2=%d\n", nos, l1, l2);
  #endif
  if (part[nos-1]) { return; }

  if (!part[l1-1]) {
    //part[l1-1]=calloc(klen+1, sizeof(double));
    computeResultRecursive(part, l1, klen);
  }
  if (!part[l2-1]) {
    //part[l2-1]=calloc(klen+1, sizeof(double));
    computeResultRecursive(part, l2, klen);
  }

#ifdef DEBUG
  Rprintf("merge l1=%d l2=%d\n", l1, l2);
  #endif
  #ifdef IN_R
  part[nos-1]=Calloc(klen+1, double);
  #else
  part[nos-1]=calloc(klen+1, sizeof(double));
  #endif
  convolute(part[nos-1], part[l1-1], part[l2-1], klen);
  #ifdef DEBUG
  for (i=0, sum=0.0;i<=klen; i++) {
    sum+=part[nos-1][i];
  }
  Rprintf("%1.3e\n ",sum);
  #endif
}

void multipleShortSequenceProbability(double *simple, double *aggregated, 
     int maxsimplehits, int maxagghits, int numofseqs) {
  int i;
  double sum;
  double **part_results;

#ifdef IN_R
  part_results=Calloc(numofseqs, double*);
  part_results[0]=Calloc(maxagghits+1, double);
  #else
  part_results=calloc(numofseqs, sizeof(double*));
  part_results[0]=calloc(maxagghits+1, sizeof(double));
  #endif
  memcpy(part_results[0],simple, (maxagghits+1)*sizeof(double));

  computeResultRecursive(part_results, numofseqs, maxagghits);

  //memcpy(aggregated, part_results[numofseqs-1], sizeof(double));

  for (i=0, sum=0.0;i<=maxagghits; i++) {
    aggregated[i]=part_results[numofseqs-1][i];
    sum+=aggregated[0];
  }
  #ifdef DEBUG
  Rprintf("P[%d]=%1.3e\n ",numofseqs, sum);
  #endif
  #ifdef IN_R
  for (i=0; i<numofseqs; i++) {
    if (part_results[i]) Free(part_results[i]);
  }
  Free(part_results);
  #else
  for (i=0; i<numofseqs; i++) {
    if (part_results[i]) free(part_results[i]);
  }
  free(part_results);
  #endif
}

