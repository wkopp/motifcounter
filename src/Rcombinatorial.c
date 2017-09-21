#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "overlap.h"
#include "combinatorial.h"
#include "score2d.h"
#include "markovchain.h"

extern double Rgran, Rsiglevel;

#define DEBUG
#undef DEBUG
void RPosteriorProbability(double *alpha, double *beta,
                           double *beta3p, double *beta5p,
                           double *hitdistribution, int *sseqlen,
                           int *smaxhits, int *snos, int *motiflen, int *singlestranded) {
    PosteriorCount prob;
    int seqlen;
    int maxhits, k, nos;
    int totalmaxhits; // total number of hits in the entire set of seq.
    double Zpartition;
    double *singlehitdistribution;
    double *delta, *deltap;
    double tau;
    double abstol = 1e-30, intol = 1e-30;
    int trace = 0, fail, fncount, type = 2, gncount;
    double sum, res;
    CGParams cgparams;


    seqlen = sseqlen[0] - motiflen[0] + 1;
    nos = snos[0];
    maxhits = smaxhits[0];
    totalmaxhits = maxhits * nos;

    delta = (double*)R_alloc((size_t)motiflen[0], sizeof(double));
    deltap = (double*)R_alloc((size_t)motiflen[0], sizeof(double));
    memset(delta, 0, motiflen[0]*sizeof(double));
    memset(deltap, 0, motiflen[0]*sizeof(double));

    computeDeltas(delta, deltap, beta, beta3p, beta5p, motiflen[0]);

    allocPosteriorProbability(&prob, seqlen, motiflen[0], maxhits);

    initPosteriorProbability(&prob, alpha[0], &beta, &beta3p, &beta5p,
                             &delta, &deltap);

    computePosteriorProbability(&prob);

    singlehitdistribution = (double*)R_alloc((size_t)(maxhits + 1),
            sizeof(double));
    memset(singlehitdistribution, 0, (maxhits + 1)*sizeof(double));

    #ifdef DEBUG
    Rprintf("omega=%e, alpha=%e\n", prob.omega, prob.alpha);
    #endif


    // renormalize the distribution so that it sums to one
    singlehitdistribution[0] = prob.probzerohits;
    Zpartition = singlehitdistribution[0];

    for (k = 1; k <= maxhits; k++) {
        finishPosteriorProbability(&prob, singlehitdistribution, k);
        Zpartition += singlehitdistribution[k];
    }
    for (k = 0; k <= maxhits; k++) {
        singlehitdistribution[k] /= Zpartition;
    }

    // compute the distribution of the number of hits
    // across multiple independent DNA sequences
    memset(hitdistribution, 0, totalmaxhits * sizeof(double));
    multipleShortSequenceProbability(singlehitdistribution, hitdistribution,
                                     maxhits, totalmaxhits, nos);
    for (k = 0, sum = 0.0; k <= totalmaxhits; k++) {
        sum += hitdistribution[k];
    }


    return;
}
