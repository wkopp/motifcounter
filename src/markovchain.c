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
#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "combinatorial.h"
#include "markovchain.h"

extern double *Rstation, *Rtrans;

double OverlapHit(int N, double *beta, double *betap) {
    int i;
    double d=1.0, n=1.0;

    //if (N<0 || N>=Rpwm->nrow) error("wrong index, i=%d\n", N);

    // beta ... forward hit
    // betap .. reverse hit either 3p or 5p
    // compute denuminator
    for (i=0; i<N; i++) {
        d-=(beta[i]+betap[i]);
    }
    n=(beta[N]);
    if (d<=0.0) return 0.0;

    return (n/d);
}

double NoOverlapHit(int N, double *beta, double *betap) {
    int i;
    double d=1.0, n=1.0;

    //if (N<0) error("wrong index, i=%d\n", N);

    // beta ... forward hit
    // betap .. reverse hit either 3p or 5p
    // compute denuminator
    for (i=0; i<N; i++) {
        d-=(beta[i]+betap[i]);
    }
    n=d-(beta[N]+betap[N]);
    if (d<=0.0) return 0.0;

    return (n/d);
}

#undef DEBUG
#define DEBUG
static double *Rdist=NULL;
void markovchain(double *dist, double *a,
 double *beta, double *beta3p, double *beta5p, int slen, int motiflen) {
    int i, k;
    double *post, *prior;
    double alphacond;

    if (!Rstation||!Rtrans) {
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
    post=Calloc(2*motiflen+2, double);
    if(post==NULL) {
        error("Memory-allocation in markovchain failed");
    }
    prior=dist;
    alphacond=a[0];
    memset(prior, 0, (2*motiflen+2)*sizeof(double));
    prior[0]=1.;

    for (k=0; k<slen; k++) {
        // P(N)
        post[0]=(1-alphacond*(2-beta3p[0]))*(prior[0]+prior[motiflen+2] + 
        prior[2*motiflen+1]);

        // P(Hf)
        post[1]=alphacond*(prior[0]+prior[motiflen+2] +
            prior[2*motiflen+1]);

        for (i=1;i<motiflen;i++) {
            post[1]+=OverlapHit(i, beta, beta3p)*prior[3+i-1];
        }
        for (i=2;i<motiflen;i++) {
            post[1]+=OverlapHit(i, beta5p, beta)*prior[motiflen+3+i-2];
        }
        post[1]+=beta5p[1]*prior[2];


        // P(Hr)
        post[2]=alphacond*(1-beta3p[0])*(prior[0]+prior[motiflen+2] +
            prior[2*motiflen+1]);
        for (i=2;i<motiflen;i++) {
            post[2]+=OverlapHit(i, beta, beta5p)*prior[motiflen+3+i-2];
        }
        for (i=1;i<motiflen;i++) {
            post[2]+=OverlapHit(i, beta3p, beta)*prior[3+i-1];
        }
        // should i switch this line
        post[2]+=beta3p[0]*prior[1];
        post[2]+=beta[1]*prior[2];

        // P(n0)
        post[3]=NoOverlapHit(0, beta, beta3p)*prior[1];
        for (i=1;i<motiflen;i++) {
            post[3+i]=NoOverlapHit(i, beta,beta3p)*prior[3+i-1];
        }
        // P(n1')
        post[3+motiflen]=NoOverlapHit(1, beta, beta5p)*prior[2];
        for (i=2;i<motiflen;i++) {
            post[motiflen+3+i-1]=NoOverlapHit(i, beta,beta5p)*
                prior[motiflen+3+i-2];
        }
        memcpy(prior, post, (2*motiflen+2)*sizeof(double));
        memset(post, 0, (2*motiflen+2)*sizeof(double));
    }

    Free(post);
}

void dmc(int n, double *alphacond, double *gradient, void *ex) {

    double val;
    CGParams *cgparams=(CGParams*)ex;
    double epsilon;
    double pa,ma;

    if (!Rstation||!Rtrans) {
        error("load forground and background properly");
        return;
    }

    if (!Rdist) {
        Rdist=Calloc(2*cgparams->motiflen+2, double);
        if(Rdist==NULL) {
            error("Memory-allocation in dmc failed");
        }
    }

    epsilon=alphacond[0]/1000;
    pa=*alphacond + epsilon;
    ma=*alphacond - epsilon;
    markovchain(Rdist, &pa, cgparams->beta, 
            cgparams->beta3p, cgparams->beta5p, 
            cgparams->len,cgparams->motiflen);

    val=Rdist[1]+Rdist[2];
    markovchain(Rdist, &ma, cgparams->beta, 
            cgparams->beta3p, cgparams->beta5p, 
            cgparams->len,cgparams->motiflen);

    val-=(Rdist[1]+Rdist[2]);
    val/=2*epsilon;

    markovchain(Rdist, alphacond, cgparams->beta, 
            cgparams->beta3p, cgparams->beta5p, 
            cgparams->len,cgparams->motiflen);

    *gradient=-2*(2*cgparams->alpha-Rdist[1]-Rdist[2])*val;
}
double minmc(int n, double *alpha, void *ex) {

    //double *extra=(double*)ex;
    CGParams *cgparams=(CGParams*)ex;

    if (!Rdist) {
        Rdist=Calloc(2*cgparams->motiflen+2, double);
        if(Rdist==NULL) {
            error("Memory-allocation in minmc failed");
        }
    }

    markovchain(Rdist, alpha, cgparams->beta, 
            cgparams->beta3p, cgparams->beta5p, 
            cgparams->len,cgparams->motiflen);

    return R_pow_di(2*cgparams->alpha-Rdist[1]-Rdist[2], 2);
}

void removeDist() {
    if(Rdist) Free(Rdist);
    Rdist=NULL;
}

