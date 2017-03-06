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
    double a0, aN;
    double abstol=1e-30, intol=1e-30;
    int trace=0, fail,fncount, type=2, gncount;
    double sum, res;
    CGParams cgparams;


    seqlen=sseqlen[0] - motiflen[0]+1;
    nos=snos[0];
    maxhits=smaxhits[0];
    totalmaxhits=maxhits*nos;

    delta=Calloc(motiflen[0],double);
    deltap=Calloc(motiflen[0],double);

    computeDeltas(delta, deltap, beta, beta3p,beta5p,motiflen[0]);

    // correct bias of alpha by fitting
    // markov model to the stationary distribution
    //extra=Calloc(3*motiflen[0]+1, double);
    //if (extra==NULL) {
        //error("Memory-allocation in RPosteriorProbability failed");
    //}
    //extra[0]=alpha[0];
    a0=alpha[0];
    //for (i=0; i<motiflen[0]; i++) {
        //extra[1+i]=beta[i];
        //extra[motiflen[0]+1+i]=beta3p[i];
        //extra[2*motiflen[0]+1+i]=beta5p[i];
    //}
    cgparams.alpha=alpha[0];
    cgparams.beta=beta;
    cgparams.beta3p=beta3p;
    cgparams.beta5p=beta5p;
    cgparams.len=500;
    cgparams.motiflen=motiflen[0];

    cgmin(1, &a0, &aN, &res, minmc, dmc, &fail, abstol, intol,
        (void*)&cgparams, type, trace, &fncount, &gncount, 100);

    //Free(extra);
    removeDist();

    allocPosteriorProbability(&prob, seqlen, motiflen[0], maxhits);
    initPosteriorProbability(&prob, alpha[0], &beta, &beta3p, &beta5p,
        &delta, &deltap);

    computePosteriorProbability(&prob);

    singlehitdistribution=Calloc(maxhits+1, double);

#ifdef DEBUG
    Rprintf("omega=%e, alpha=%e\n", prob.omega, prob.alpha);
#endif


    // renormalize the distribution so that it sums to one
    singlehitdistribution[0]=prob.probzerohits;
    Zpartition=singlehitdistribution[0];

    for (k=1; k<=maxhits; k++) {
        finishPosteriorProbability(&prob, singlehitdistribution, k);
        Zpartition+=singlehitdistribution[k];
    }
    for (k=0; k<=maxhits; k++) {
        singlehitdistribution[k]/=Zpartition;
    }

    // compute the distribution of the number of hits
    // across multiple independent DNA sequences
    multipleShortSequenceProbability(singlehitdistribution, hitdistribution,
                maxhits, totalmaxhits, nos);
    for (k=0,sum=0.0; k<=totalmaxhits; k++) {
        sum+=hitdistribution[k];
    }

    deletePosteriorProbability(&prob);
    Free(singlehitdistribution);
    Free(delta);
    Free(deltap);

    return;
}
