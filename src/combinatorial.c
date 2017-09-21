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
#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
//#include "countdist.h"
#include "combinatorial.h"
#include "markovchain.h"


int allocPosteriorProbability(PosteriorCount *p, int seqlen,
                              int mlen, int maxhits) {
    int i, j;
    p->seqlen = seqlen;
    p->mlen = mlen;
    p->maxhits = maxhits;
    p->value = (double***)R_alloc((size_t)maxhits, sizeof(double **));
    memset(p->value, 0, maxhits * sizeof(double **));

    for (i = 0; i < maxhits; i++) {
        p->value[i] = (double**)R_alloc((size_t)seqlen, sizeof(double *));
        memset(p->value[i], 0, seqlen * sizeof(double *));
        for (j = 0; j < seqlen; j++) {
            p->value[i][j] = (double*)R_alloc(2 * mlen, sizeof(double));
            memset(p->value[i][j], 0, 2 * mlen * sizeof(double));
        }
    }
    return 0;
}

double addomegas(double *omega, int start, int end) {
    int i;
    double x = 1.0;

    for (i = start; i <= end; i++) {
        x *= omega[i];
    }
    return x;
}

void printVector(double **m, int i1, int i2, int len) {
#ifndef IN_R
    int i;
    for (i = 0; i < len; i++) {
        printf("%1.2e\t", m[i1][i2 + i]);
    }
    printf("\n");
#endif
}

#define DEBUG
#undef DEBUG
void initPosteriorProbability(PosteriorCount *p, double alpha, double **beta,
                              double **beta3p, double **beta5p, double **delta, double **deltap) {
    //double *extra;
    int i, j;
    int m;
    double abstol = 1e-30, intol = 1e-30;
    double res;
    int trace = 0, fail, fncount, type = 2, gncount;
    CGParams cgparams;
    double a0, aN;
    double *_alpha, *_omega;
    double tau;

    p->beta = *beta;
    p->beta3p = *beta3p;
    p->beta5p = *beta5p;
    p->delta = *delta;
    p->deltap = *deltap;

    _alpha = (double*)R_alloc((size_t)p->seqlen, sizeof(double));
    _omega = (double*)R_alloc((size_t)p->seqlen, sizeof(double));
    memset(_alpha, 0, (size_t)p->seqlen*sizeof(double));
    memset(_omega, 0, (size_t)p->seqlen*sizeof(double));

    m = (70 > p->seqlen) ? p->seqlen : 70;


    #ifdef WKO
    a0 = alpha;

    cgparams.alpha = alpha;
    cgparams.beta = p->beta;
    cgparams.beta3p = p->beta3p;
    cgparams.beta5p = p->beta5p;
    cgparams.len = 500;
    cgparams.motiflen = p->mlen;
    cgparams.dist = (double*)R_alloc((size_t)2 * cgparams.motiflen + 2,
            sizeof(double));
    memset(cgparams.dist, 0, (2 * cgparams.motiflen+2)*sizeof(double));
    /*
       for (i=0; i<p->mlen; i++) {
       extra[1+i]=p->beta[i];
       extra[p->mlen+1+i]=p->beta3p[i];
       extra[2*p->mlen+1+i]=p->beta5p[i];
       }
       */

    m = 1;
    for (i = 0; i < m; i++) {

        //   extra[3*p->mlen+1]=(double)(500);

        cgmin(1, &a0, &aN, &res, minmc, dmc, &fail, abstol, intol,
              (void *)&cgparams, type, trace, &fncount, &gncount, 100);

        _alpha[i] = tau;
        a0 = aN;
    }
    #else

    tau = getOptimalTauMCDS(&alpha, p->beta, p->beta3p, p->beta5p, &p->mlen);

    m = 1;
    for (i = 0; i < m; i++) {
        _alpha[i] = tau;
    }

    #endif

    for (i = m; i < p->seqlen; i++) {
        _alpha[i] = _alpha[m - 1];
    }
    for (i = 0; i < p->seqlen; i++) {
        _omega[i] = 1 - 2 * _alpha[i] + _alpha[i] * p->beta3p[0];
    }
    p->probzerohits = addomegas(_omega, 0, p->seqlen - 1);


    p->alpha = tau;
    p->omega = 1 - 2 * p->alpha + p->alpha * p->beta3p[0];

#ifdef DEBUG
#ifdef IN_R
    Rprintf("a=%e, b=%e,b3p=%e, b5p=%e, d=%e, dp=%e, o=%e\n",
            alpha, p->beta[1], p->beta3p[0], p->beta5p[1], p->delta[0], p->deltap[0],
            p->omega);
#endif
#endif

    // init k=1
    for (j = 0; j < p->seqlen; j++) {
        for (i = 0; i <= j; i++) {
            if (p->mlen - 1 <= j - i) {

                p->value[0][j][0] +=
                                _alpha[i] * p->delta[p->mlen - 1]
                                * addomegas(_omega, 0, i - 1)
                                * addomegas(_omega, i + p->mlen, j);

                p->value[0][j][p->mlen] +=
                                _alpha[i] * p->delta[0] * p->deltap[p->mlen - 1]
                                * addomegas(_omega, 0, i - 1)
                                * addomegas(_omega, i + p->mlen, j);

            } else {
                p->value[0][j][p->mlen - j + i - 1] +=
                                _alpha[i] * addomegas(_omega, 0, i - 1);
                p->value[0][j][2 * p->mlen - j + i - 1] += _alpha[i] *
                        p->delta[0] * addomegas(_omega, 0, i - 1);
            }
        }

#ifdef DEBUG
        for (i = 0; i < 2 * p->mlen; i++) {
            Rprintf( "%1.2e\t", p->value[0][j][i]);
        }
        Rprintf( "\n");
#endif
    }
}

double ffTransProb(PosteriorCount *prob, int n) {
    if (n < 0 ) return 0.0;
    else if (n < prob->mlen) return (prob->beta[n]);
    else return prob->alpha;
}

double rrTransProb(PosteriorCount *prob, int n) {
    if (n < 0 ) return 0.0;
    else if (n < prob->mlen) return (prob->beta[n]);
    else return prob->alpha * prob->delta[0];
}

double frTransProb(PosteriorCount *prob, int n) {
    if (n < 0 ) return 0.0;
    else if (n < prob->mlen) return (prob->beta3p[n]);
    else return prob->alpha * prob->delta[0];
}

double rfTransProb(PosteriorCount *prob, int n) {
    if (n < 0 ) return 0.0;
    else if (n < prob->mlen) return (prob->beta5p[n]);
    else return prob->alpha ;
}

double rNonHitStretch(PosteriorCount *prob, int n) {
    if (n < prob->mlen) return 1.0;
    else {
        return prob->deltap[prob->mlen - 1] * R_pow_di(prob->omega, n - prob->mlen);
    }
}

double fNonHitStretch(PosteriorCount *prob, int n) {
    if (n < prob->mlen) return 1.0;
    else return prob->delta[prob->mlen - 1] * R_pow_di(prob->omega,
                    n - prob->mlen);
}

#define DEBUG
#undef DEBUG
void computePosteriorProbability(PosteriorCount *prob) {
    int i, j, k;
    int mc, mp;
#ifdef DEBUG
    int p;
#endif
    for (k = 1; k < prob->maxhits; k++) {

#ifndef IN_R
#ifdef DEBUG
        fprintf(stdout, "k=%d ... \n", k);
#endif
#endif

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(k,prob) private(i,j, mc,mp)
#endif
        for (i = 0; i < prob->seqlen; i++) {
            R_CheckUserInterrupt();
            for (j = 0; j < i; j++) {

                if (prob->mlen <= i - j) mc = 0;
                else mc = prob->mlen - i + j;

                for (mp = 0; mp < prob->mlen; mp++) {
                    // forward hit
                    prob->value[k][i][mc] +=
                                    (prob->value[k - 1][j][mp] *
                                     ffTransProb(prob, prob->mlen - mp) +
                                     prob->value[k - 1][j][prob->mlen + mp] *
                                     rfTransProb(prob, prob->mlen - mp))
                                    * fNonHitStretch(prob, i - j);

                    // reverse hit
                    prob->value[k][i][prob->mlen + mc] +=
                                    (prob->value[k - 1][j][prob->mlen + mp] *
                                     rrTransProb(prob, prob->mlen - mp) +
                                     prob->value[k - 1][j][mp] *
                                     frTransProb(prob, prob->mlen - mp))
                                    * rNonHitStretch(prob, i - j);


                    ////////////////////////
                    //
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
                    fprintf(stdout, "p(X_[0,%d],a=%d|k=%d)+=(p([0,%d],%d|%d)*"
                            "%1.1e + p([0,%d]',%d|%d)*%1.1e)*%1.1e=\n ",
                            i, mc, k,
                            j, mp, k - 1,
                            ffTransProb(prob, prob->mlen - mp),
                            j, mp, k - 1,
                            rfTransProb(prob, prob->mlen - mp), fNonHitStretch(prob, i - j));

                    fprintf(stdout, "\t (%1.1e*%1.1e+%1.1e*%1.1e)*%1.1e\n",
                            prob->value[k - 1][j][mp],
                            ffTransProb(prob, prob->mlen - mp),
                            prob->value[k - 1][j][prob->mlen + mp],
                            rfTransProb(prob, prob->mlen - mp), fNonHitStretch(prob, i - j));

                    fprintf(stdout, "p(X_[0,%d]',a=%d|k=%d)+=(p([0,%d]',%d|%d)*"
                            "%1.1e + p([0,%d],%d|%d)*%1.1e)*%1.1e=\n ",
                            i, mc, k,
                            j, mp, k - 1,
                            rrTransProb(prob, prob->mlen - mp),
                            j, mp, k - 1,
                            frTransProb(prob, prob->mlen - mp), fNonHitStretch(prob, i - j));

                    fprintf(stdout, "\t (%1.1e*%1.1e+%1.1e*%1.1e)*%1.1e\n",
                            prob->value[k - 1][j][prob->mlen + mp],
                            rrTransProb(prob, prob->mlen - mp),
                            prob->value[k - 1][j][mp],
                            frTransProb(prob, prob->mlen - mp), rNonHitStretch(prob, i - j));
#endif
#endif

                }
            }

            for (j = 0; j <= i; j++) {
                R_CheckUserInterrupt();
                if (prob->mlen - 1 <= i - j) mc = 0;
                else mc = prob->mlen - i + j - 1;

                prob->value[k][i][prob->mlen + mc] +=
                                prob->value[k - 1][j][prob->mlen - 1] *
                                frTransProb(prob, 0) * rNonHitStretch(prob, i - j + 1);

#define DEBUG
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
                fprintf(stdout, "p(X_[0,%d]',a=%d|k=%d)+="
                        "p([0,%d],%d|%d) * %1.1e * %1.1e=\n ",
                        i, mc, k,
                        j, prob->mlen - 1, k - 1,
                        frTransProb(prob, 0), rNonHitStretch(prob, i - j + 1));

                fprintf(stdout, "\t %1.1e * %1.1e * %1.1e\n",
                        prob->value[k - 1][j][prob->mlen - 1],
                        frTransProb(prob, 0), rNonHitStretch(prob, i - j + 1));
#endif
#endif
            }

#define DEBUG
#undef DEBUG
#ifndef IN_R
#ifdef DEBUG
            for(p = 0; p < prob->mlen * 2; p++) {
                fprintf(stdout, "%1.2e\t", prob->value[k][i][p]);
            }
            fprintf(stdout, "\n");
#endif
#endif
        }

    }
}

void finishPosteriorProbability(PosteriorCount *prob,
                                double *final, int nhits) {
    int m;
    final[nhits] += prob->value[nhits - 1][prob->seqlen - 1][0];
    for (m = 1; m < prob->mlen; m++) {
        final[nhits] += prob->value[nhits - 1][prob->seqlen - 1][m] *
                        (prob->delta[prob->mlen - m - 1]);
    }
    final[nhits] += prob->value[nhits - 1][prob->seqlen - 1][prob->mlen];
    for (m = 1; m < prob->mlen; m++) {
        final[nhits] += prob->value[nhits - 1][prob->seqlen - 1][prob->mlen + m] *
                        (prob->deltap[prob->mlen - m - 1]);
    }
}

#undef DEBUG
// convolution operation
void convolute(double *result, double *p1, double *p2, int len) {
    int i, j;
    for(i = 0; i <= len; i++) {
        for (j = 0; j <= len && i + j <= len; j++) {
            result[i + j] += p1[i] * p2[j];
        }
    }
}

// this function determines the actual recursive convolution
// of the individual sequences.
// The final result is determined from intermediate results
// in a dynamic programming manner to reduce redundant computational
// effort.
void computeResultRecursive(double **part, int nos, int klen) {
    int l1, l2;
#ifdef DEBUG
    int i;
    double sum;
#endif
    l1 = nos / 2;
    l2 = nos - l1;

#ifdef DEBUG
    Rprintf("nos=%d from l1=%d l2=%d\n", nos, l1, l2);
#endif
    // if the distribution across nos sequences was reached,
    // the computation is accomplished, otherwise continue
    // the recursion
    if (part[nos - 1]) {
        return;
    }

    // if the first half of the computation wasn't
    // evaluated previously, call computeResultRecursive
    // with l1
    if (!part[l1 - 1]) {
        computeResultRecursive(part, l1, klen);
    }
    // if the second half of the computation wasn't
    // evaluated previously, call computeResultRecursive
    // with l2
    if (!part[l2 - 1]) {
        computeResultRecursive(part, l2, klen);
    }

#ifdef DEBUG
    Rprintf("merge l1=%d l2=%d\n", l1, l2);
#endif
    part[nos - 1] = (double*)R_alloc((size_t)klen + 1, sizeof(double));
    memset(part[nos - 1], 0, (klen + 1)*sizeof(double));
    convolute(part[nos - 1], part[l1 - 1], part[l2 - 1], klen);
#ifdef DEBUG
    for (i = 0, sum = 0.0; i <= klen; i++) {
        sum += part[nos - 1][i];
    }
    Rprintf("%1.3e\n ", sum);
#endif
}

// This function determines the distribution of the number of hits
// in multiple i.i.d. sequences based on recursively
// convolving the distribution of hits of single sequences
void multipleShortSequenceProbability(double *singledist,
                                      double *aggregateddist,
                                      int maxsinglehits, int maxagghits, int numofseqs) {
    int i;
    double **part_results;

    // allocate array of pointers to sub-aggregated distributions
    part_results = (double**)R_alloc((size_t)numofseqs, sizeof(double *));
    memset(part_results, 0, numofseqs * sizeof(double *));

    part_results[0] = (double*)R_alloc((size_t)maxagghits + 1, sizeof(double));

    //memset(part_results[0], 0, (maxagghits + 1) * sizeof(double));

    // copy the distribution of a single sequence to part_results
    // this corresponds to the base case
    memset(part_results[0], 0, (maxagghits+1) * sizeof(double));
    memcpy(part_results[0], singledist, (maxsinglehits + 1)*sizeof(double));

    // recursively determine the distribution in N i.i.d. sequences
    computeResultRecursive(part_results, numofseqs, maxagghits);

    // copy the final result into aggregated
    for (i = 0; i <= maxagghits; i++) {
        aggregateddist[i] = part_results[numofseqs - 1][i];
    }
}
