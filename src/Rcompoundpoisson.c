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

    seqlen = 0;
    for (i = 0; i < *nseq; i++) {
        seqlen += lseq[i] - motiflen[0] + 1;
    }

    maxclumpsize = (double)mclump[0];
    maxhits = (double)mhit[0];
    singlestranded = *sstrand;

    if (singlestranded == 1) {

        theta = initThetaSingleStranded(maxclumpsize);

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
    double *theta;
    int i;
  
    theta = initTheta(maxclump[0]);
  
    clumpsizeBeta(beta, beta3p, beta5p, theta, maxclump, motiflen);
    
    for (i = 0; i < maxclump[0]; i++) {
      dist[i] = theta[i*2] + theta[i*2 + 1];
    }
}

void RclumpsizeGamma(double *gamma, double *dist, int *maxclump, int *motiflen) {
    double *theta;
    int i;
    
    theta = initTheta(maxclump[0]);
  
    clumpsizeGamma(gamma, theta, maxclump, motiflen);
    
    for (i = 0; i < maxclump[0]; i++) {
      dist[i] = theta[i*2] + theta[i*2 + 1];
    }
}

