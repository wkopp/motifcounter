#ifdef IN_R
#include <R.h>
#endif

#include "matrix.h"
#include "score2d.h"
#include "overlap.h"
#include "compoundpoisson.h"

#define EPSILON 1e-10
extern double Rsiglevel, Rgran;

void RcompoundpoissonPape_useGamma(double *gamma,
    double *hitdistribution, int *nseq, int *lseq, int * mhit, int *mclump,
    int *motiflen) {
    int seqlen, i;
    int maxclumpsize, maxhits;
    double lambda;
    double *theta, extention[3];

    seqlen=0;
    for (i=0; i<*nseq; i++) {
        if (lseq[i]-motiflen[0]+1 > 0) {
            seqlen+=lseq[i]-motiflen[0]+1;
        }
    }
    maxclumpsize=(double)mclump[0];
    maxhits=(double)mhit[0];

    memset(extention, 0, 3*sizeof(double));
    gamma[motiflen[0]]=(gamma[motiflen[0]]+EPSILON)/(1+2*EPSILON);
    computeExtentionFactorsPape(extention, gamma, motiflen[0]);
    theta=initTheta(maxclumpsize);

    computeInitialClump(theta, gamma,motiflen[0]);
    computeTheta(maxclumpsize, theta, extention, motiflen[0]);

    lambda=computePoissonParameter(seqlen, motiflen[0],
            maxclumpsize, gamma[0],theta);
    computeCompoundPoissonDistributionKemp(lambda, maxhits,
            maxclumpsize, theta, hitdistribution);
    deleteTheta(theta);
}

void Rcompoundpoisson_useBeta(double *alpha, double *beta,
                double *beta3p, double *beta5p,
                double *hitdistribution, int *nseq, 
                int *lseq, int * mhit, int *mclump, int *motiflen,
                int *sstrand) {
    int seqlen, i;
    int maxclumpsize, maxhits, singlestranded;
    double lambda;
    double *theta, extention[3];
    double *delta, *deltap;

    //compute the total length of the sequence
    seqlen=0;
    for (i=0; i<*nseq; i++) {
        seqlen+=lseq[i]-motiflen[0]+1;
    }
    //Rprintf("nseq=%d, seqlen=%d\n",*nseq,seqlen);
    //init the maximal clump size and the max number of hits
    maxclumpsize=(double)mclump[0];
    maxhits=(double)mhit[0];
    singlestranded=*sstrand;

    delta=Calloc(motiflen[0],double);
    deltap=Calloc(motiflen[0],double);

    // initialize the extention factors
    memset(extention, 0, 3*sizeof(double));

    if (singlestranded==1) {
        computeDeltasSingleStranded(delta, beta, motiflen[0]);

        computeExtentionFactorsKoppSingleStranded(extention, beta, motiflen[0]);
        theta=initThetaSingleStranded(maxclumpsize);

        computeInitialClumpKoppSingleStranded(theta, delta, motiflen[0]);
        computeThetaSingleStranded(maxclumpsize, theta, extention, motiflen[0]);

        lambda=computePoissonParameterSingleStranded(seqlen, motiflen[0],
                    maxclumpsize, alpha[0],theta);

        computeCompoundPoissonDistributionKempSingleStranded(lambda, maxhits,
                    maxclumpsize, theta, hitdistribution);
    } else {
        beta3p[0]=(beta3p[0]+EPSILON)/(1+2*EPSILON);
        computeDeltas(delta, deltap, beta, beta3p,beta5p,motiflen[0]);


        computeExtentionFactorsKopp(extention, delta, deltap, beta,
                    beta3p, beta5p, motiflen[0]);
        theta=initTheta(maxclumpsize);

        computeInitialClumpKopp(theta, beta3p,delta, deltap, motiflen[0]);
        computeTheta(maxclumpsize, theta, extention, motiflen[0]);

        lambda=computePoissonParameter(seqlen, motiflen[0], maxclumpsize,
                                alpha[0],theta);
        //Rprintf("lambda=%e\n",lambda);
        computeCompoundPoissonDistributionKemp(lambda, maxhits, maxclumpsize,
                                theta, hitdistribution);
    }



    deleteTheta(theta);
    Free(delta);
    Free(deltap);
}
