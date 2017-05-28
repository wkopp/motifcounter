#ifdef IN_R
#include <R.h>
#endif

#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "compoundpoisson.h"

extern double Rsiglevel, Rgran;

void RcompoundpoissonPape_useGamma(double *gamma,
                                   double *hitdistribution, int *nseq, int *lseq, int *mhit, int *mclump,
                                   int *motiflen) {
    int seqlen, i;
    int maxclumpsize, maxhits;
    double lambda;
    double *theta, extention[3];

    seqlen = 0;
    for (i = 0; i < *nseq; i++) {
        if (lseq[i] - motiflen[0] + 1 > 0) {
            seqlen += lseq[i] - motiflen[0] + 1;
        }
    }
    maxclumpsize = (double)mclump[0];
    maxhits = (double)mhit[0];
    theta = initTheta(maxclumpsize);

    clumpsizeGamma(gamma, theta, mclump, motiflen);

    lambda = computePoissonParameter(seqlen, motiflen[0],
                                     maxclumpsize, gamma[0], theta);
    computeCompoundPoissonDistributionKemp(lambda, maxhits,
                                           maxclumpsize, theta, hitdistribution);
}

void Rcompoundpoisson_useBeta(double *alpha, double *beta,
                              double *beta3p, double *beta5p,
                              double *hitdistribution, int *nseq,
                              int *lseq, int *mhit, int *mclump, int *motiflen,
                              int *sstrand) {
    int seqlen, i;
    int maxclumpsize, maxhits, singlestranded;
    double lambda;
    double *theta, extention[3];
    double *delta, *deltap;

    //compute the total length of the sequence
    seqlen = 0;
    for (i = 0; i < *nseq; i++) {
        seqlen += lseq[i] - motiflen[0] + 1;
    }
    //Rprintf("nseq=%d, seqlen=%d\n",*nseq,seqlen);
    //init the maximal clump size and the max number of hits
    maxclumpsize = (double)mclump[0];
    maxhits = (double)mhit[0];
    singlestranded = *sstrand;

    //delta = (double*)R_alloc((size_t)motiflen[0], sizeof(double));
    //deltap = (double*)R_alloc((size_t)motiflen[0], sizeof(double));

    //memset(delta, 0, motiflen[0]*sizeof(double));
    //memset(deltap, 0, motiflen[0]*sizeof(double));

    // initialize the extention factors
    //memset(extention, 0, 3 * sizeof(double));

    if (singlestranded == 1) {
        //computeDeltasSingleStranded(delta, beta, motiflen[0]);

        //computeExtentionFactorsKoppSingleStranded(extention, beta, motiflen[0]);
        theta = initThetaSingleStranded(maxclumpsize);

        //computeInitialClumpKoppSingleStranded(theta, delta, motiflen[0]);
        //computeThetaSingleStranded(maxclumpsize, theta, extention, motiflen[0]);
        clumpsizeBeta_singlestranded(beta, theta, mclump, motiflen);

        lambda = computePoissonParameterSingleStranded(seqlen, motiflen[0],
                 maxclumpsize, alpha[0], theta);

        computeCompoundPoissonDistributionKempSingleStranded(lambda, maxhits,
                maxclumpsize, theta, hitdistribution);
    } else {
        theta = initTheta(maxclumpsize);
        clumpsizeBeta(beta, beta3p, beta5p, theta, mclump, motiflen);

        lambda = computePoissonParameter(seqlen, motiflen[0], maxclumpsize,
                                         alpha[0], theta);
        computeCompoundPoissonDistributionKemp(lambda, maxhits, maxclumpsize,
                                               theta, hitdistribution);
    }
}

void RclumpsizeBeta(double *beta, double *beta3p, double *beta5p,
                        double *dist, int *maxclump, int *motiflen) {
    clumpsizeBeta(beta, beta3p, beta5p, dist, maxclump, motiflen);
}

void RclumpsizeGamma(double *gamma, double *dist, int *maxclump, int *motiflen) {

    clumpsizeGamma(gamma, dist, maxclump, motiflen);
}

