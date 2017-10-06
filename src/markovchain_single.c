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
#endif

#include "background.h"
#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "combinatorial.h"
#include "markovchain.h"


static double OverlapHit(int N, double *beta) {
    int i;
    double d = 1.0, n = 1.0;

    //if (N<0 || N>=Rpwm->nrow) error("wrong index, i=%d\n", N);

    // beta ... forward hit
    // betap .. reverse hit either 3p or 5p
    // compute denuminator
    for (i = 0; i < N; i++) {
        d -= (beta[i]);
    }
    n = (beta[N]);
    if (d <= 0.0) return 0.0;

    return (n / d);
}

static double NoOverlapHit(int N, double *beta) {
    int i;
    double d = 1.0, n = 1.0;

    //if (N<0) error("wrong index, i=%d\n", N);

    // beta ... forward hit
    // betap .. reverse hit either 3p or 5p
    // compute denuminator
    for (i = 0; i < N; i++) {
        d -= (beta[i]);
    }
    n = d - (beta[N]);
    if (d <= 0.0) return 0.0;

    return (n / d);
}

#undef DEBUG
#define DEBUG
void markovchain_ss(double *dist, double *tau_,
                 double *beta, int *slen_, int *motiflen_) {
    int i, k;
    int slen = *slen_;
    int motiflen = *motiflen_;
    double *post, *prior;
    double tau = tau_[0];


    // the states are
    // dist[0] ... p(nohit)
    // dist[1] ... p(H)
    // dist[2 ... M-2] ... p(n0), ... , p(nL)
    //
    //Rprintf("motiflen=%d, slen=%d", motiflen, slen);

    post = (double*)R_alloc((size_t)motiflen, sizeof(double));
    memset(post, 0, (motiflen)*sizeof(double));

    prior = dist;
    memset(prior, 0, (motiflen)*sizeof(double));
    prior[0] = 1.;

    for (k = 0; k < slen; k++) {
        // P(N)
        post[0] = (1 - tau) * (prior[0] + prior[motiflen - 1]);

        // P(Hf)
        post[1] = tau * (prior[0] + prior[motiflen - 1]);

        for (i = 0; i < motiflen - 2; i++) {
            post[1] += OverlapHit(i + 1, beta) * prior[1 + i];
        }

        // P(n0)
        //post[2] = NoOverlapHit(0, beta) * prior[1];
        for (i = 0; i < motiflen - 2; i++) {
            post[2 + i] = NoOverlapHit(i + 1, beta) * prior[2 + i - 1];
        }

        memcpy(prior, post, (motiflen)*sizeof(double));
        memset(post, 0, (motiflen)*sizeof(double));
    }
}

static double minmc_ss(int n, double *tau, void *ex) {

    CGParams *cgparams = (CGParams *)ex;

    markovchain_ss(cgparams->dist, tau, cgparams->beta,
                &cgparams->len, &cgparams->motiflen);

    return - cgparams->alpha * log(cgparams->dist[1]) -
            (1 - cgparams->alpha) * log(1 - cgparams->dist[1]);
}


static void dmc_ss(int n, double *tau, double *gradient, void *ex) {

    double val;
    double epsilon;
    double pa, ma;

    epsilon = tau[0] / 1000;
    pa = *tau + epsilon;
    ma = *tau - epsilon;

    val = (minmc_ss(n, &pa, ex) - minmc_ss(n, &ma, ex)) / (2 * epsilon);

    *gradient =  val;
}


// Returns the clump start probabilíty for the given markov model
double getOptimalTauMCSS(double *alpha, double *beta, int *motiflen) {
    double a0, aN;
    double abstol = 1e-30, intol = 1e-30;
    int trace = 0, fail, fncount, type = 2, gncount;
    double res;
    CGParams cgparams;

    a0 = alpha[0];
    cgparams.alpha = alpha[0];
    cgparams.beta = beta;
    cgparams.len = 500;
    cgparams.motiflen = motiflen[0];
    cgparams.dist = (double*)R_alloc((size_t) cgparams.motiflen,
            sizeof(double));


    cgmin(1, &a0, &aN, &res, minmc_ss, dmc_ss, &fail, abstol, intol,
          (void *)&cgparams, type, trace, &fncount, &gncount, 100);

    return a0;
}


SEXP mcss_check_optimal(SEXP alpha_, SEXP beta_, SEXP motiflen_) {

    double *alpha = REAL(alpha_);
    double *beta = REAL(beta_);
    int *motiflen = INTEGER(motiflen_);

    double tau;

    tau = getOptimalTauMCSS(alpha, beta, motiflen);

    return ScalarReal(tau);
}
