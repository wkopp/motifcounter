#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef IN_R
#include <R.h>
#include <Rmath.h>
#endif
#include "score1d.h"
#include "score2d.h"
#include "compoundpoisson.h"
#include "background.h"
#include "overlap.h"

#define DSTRANDED 2
#define EPSILON 1e-10
// This function determines the Poisson parameter
// for scanning both strands of the DNA sequence
//
//      2 x alpha x (S-M+1)
//      -------------------
//         E[ clumpsize]
//
double computePoissonParameter(int seqlen, int mlen,
                               int maxclump, double alpha, double *theta) {
    int i;
    double ec = 0.0;
    for (i = 0; i < maxclump; i++) {
        ec += (theta[i * DSTRANDED] + theta[i * DSTRANDED + 1]) * (double)(i + 1);
    }
    return ((double)2 * (seqlen) * (alpha)) / (ec);
}

// This function determines the Poisson parameter
// for scanning a single strand of the DNA sequence
//
//        alpha x (S-M+1)
//      -------------------
//         E[ clumpsize]
//
double computePoissonParameterSingleStranded(int seqlen, int mlen,
        int maxclump, double alpha, double *theta) {
    int i;
    double ec = 0.0;
    for (i = 0; i < maxclump; i++) {
        ec += (theta[i]) * (double)(i + 1);
    }
    return ((double)(seqlen) * (alpha)) / (ec);
}

// allocate the memory required to compute the clump size distribution
double *initTheta(int maxclump) {
    double *p;
    p = (double*)R_alloc((size_t)maxclump * DSTRANDED, sizeof(double));
    memset(p, 0, (size_t)maxclump * DSTRANDED* sizeof(double));
    return p;
}

double *initThetaSingleStranded(int maxclump) {
    double *p;

    p = (double*)R_alloc((size_t)maxclump, sizeof(double));
    memset(p, 0, (size_t)maxclump * sizeof(double));
    return p;
}

// This is the Pape et al. version of having a 1-clump
void computeInitialClump(double *theta, double *gamma, int mlen) {
    int i;
    theta[0] = 1 - gamma[mlen];
    theta[1] = 1 - gamma[mlen];

    for (i = 1; i < mlen; i++) {
        theta[0] *= (1 - gamma[i]) * (1 - gamma[mlen + i]);
        theta[1] *= (1 - gamma[i]) * (1 - gamma[mlen * 2 + i]);
    }
}

// This is the Kopp et al. version of having a 1-clump for
// scanning both strands of the DNA sequence
void computeInitialClumpKopp(double *theta, double *beta3p,
                             double *delta, double *deltap, int mlen) {
    theta[0] = delta[mlen - 1];
    theta[1] = deltap[mlen - 1] * (1 - beta3p[0]);
}

// This is the Kopp et al. version of having a 1-clump for
// scanning a single strand of the DNA sequence
void computeInitialClumpKoppSingleStranded(double *theta,
        double *delta,  int mlen) {
    theta[0] = delta[mlen - 1];
}

// Computation of the extention factors according to Pape et al.
// Note that this function is only applicable to scanning both DNA strands
void computeExtentionFactorsPape(double *xi, double *gamma, int mlen) {
    int k, j;
    double xik;
    // xi
    for (k = 1; k < mlen; k++) {
        xik = gamma[k] / (1 - gamma[k]);
        xik *= (1 - gamma[mlen]) / (1 - gamma[mlen + k]);
        for (j = 1; j < mlen - k; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen + j]) /
                   ((1 - gamma[k + j]) * (1 - gamma[mlen + k + j]));
        }
        for (j = mlen - k; j < mlen; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen + j]);
        }
#ifdef DEBUG
        Rprintf( "xi%d=%f\n", k, xik);
#endif
        xi[0] += xik;
    }

    // xi',3'
    for (k = 0; k < mlen; k++) {
        if (k == 0) {
            xik = (gamma[mlen + k]) / (1 - gamma[mlen + k]);
        } else {
            xik = gamma[mlen + k] / (1 - gamma[mlen + k]);
        }
        for (j = 1; j < mlen - k; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen * 2 + j]) /
                   ((1 - gamma[k + j]) * (1 - gamma[mlen + k + j]));
        }
        for (j = mlen - k; j < mlen; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen * 2 + j]);
        }
#ifdef DEBUG
        Rprintf( "xi%d'=%f\n", k, xik);
#endif
        xi[1] += xik;
    }

    // xi,5'
    for (k = 1; k < mlen; k++) {
        xik = gamma[mlen * 2 + k] / (1 - gamma[mlen * 2 + k]);
        xik *= (1 - gamma[mlen]) / (1 - gamma[k]);
        for (j = 1; j < mlen - k; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen + j]) /
                   ((1 - gamma[k + j]) * (1 - gamma[mlen * 2 + k + j]));
        }
        for (j = mlen - k; j < mlen; j++) {
            xik *= (1 - gamma[j]) * (1 - gamma[mlen + j]);
        }
#ifdef DEBUG
        Rprintf ("xi%d,5'=%f\n", k, xik);
#endif
        xi[2] += xik;
    }
#ifdef DEBUG
    Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n", xi[0], xi[1], xi[2]);
#endif
}

// This function computes the extension factors according to Kopp et al.
// for scanning both DNA strands
void computeExtentionFactorsKopp(double *xi,
                                 double *delta, double *deltap, double *beta,
                                 double *beta3p, double *beta5p, int mlen) {
    int k;

    xi[1] = beta3p[0];
    for (k = 1; k < mlen; k++) {
        xi[0] += beta[k];
        xi[1] += beta3p[k];
        xi[2] += beta5p[k];
    }
    xi[1] *= deltap[mlen - 1] / (delta[mlen - 1]);
    xi[2] *= (delta[mlen - 1]) / deltap[mlen - 1];
#ifdef DEBUG
    Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n", xi[0], xi[1], xi[2]);
#endif

}

// This function computes the extension factors according to Kopp et al.
// for scanning a single DNA strand
void computeExtentionFactorsKoppSingleStranded(double *xi, double *beta,
        int mlen) {
    int k;

    for (k = 1; k < mlen; k++) {
        xi[0] += beta[k];
    }
#ifdef DEBUG
    Rprintf("Extention factors: xi=%e\n", xi[0]);
#endif

}

void computeTheta(int maxclump, double *theta,
                  double *extention, int mlen) {
    int i;
    double total = 0.0;

    total = theta[0] + theta[1];
    for (i = 1; i < maxclump; i++) {
        theta[i * DSTRANDED] =
                        extention[0] * theta[(i - 1) * DSTRANDED] +
                        extention[2] * theta[(i - 1) * DSTRANDED + 1];
        theta[i * DSTRANDED + 1] =
                        extention[1] * theta[(i - 1) * DSTRANDED] +
                        extention[0] * theta[(i - 1) * DSTRANDED + 1];

        total += (theta[i * DSTRANDED] + theta[i * DSTRANDED + 1]);

    }

    for (i = 0; i < maxclump; i++) {
        theta[i * DSTRANDED] /= total;
        theta[i * DSTRANDED + 1] /= total;
    }
}

void computeThetaSingleStranded(int maxclump, double *theta,
                                double *extention, int mlen) {
    int i;
    double total = 0.0;

    // geometric series appoximation
    total = theta[0];
    for (i = 1; i < maxclump; i++) {
        theta[i] = extention[0] * theta[(i - 1)];

        total += (theta[i]);
    }

    // normalize the theta such that they sum up to one
    for (i = 0; i < maxclump; i++) {
        theta[i] /= total;
    }
}


// Implementation of Kemp et al., which was also used
// in Pape et al.
void computeCompoundPoissonDistributionKemp(double lambda,
        int maxhit, int maxclump, double *theta, double *cp) {
    int i, j;
    double p;
    double normalize = 0.0;
    double logminp, maxcp;

#ifdef DEBUG
    Rprintf( "lambda=%f\n", lambda);
#endif
    // start log(cp[0]))
    cp[0] = -lambda;

    // compute the compound poisson distribution according to
    // Kemp in log-space
    for (i = 1; i <= maxhit; i++) {
        p = 0.0;
        j = i - maxclump + 1;
        j = (j > 0) ? j : 0;
        // find the smallest logp
        logminp = cp[j];
        for (; j < i; j++) {
            if (logminp > cp[j]) {
                logminp = cp[j];
            }
        }
        j = i - maxclump + 1;
        j = (j > 0) ? j : 0;
        // compute p which is divided by logminp
        for (; j < i; j++) {
            p += (i - j) * (theta[(i - j - 1) * DSTRANDED] + theta[(i - j - 1) * DSTRANDED
                            + 1]) *
                 exp(cp[j] - logminp);
        }
        // compute log(cp)
        cp[i] = log(lambda / ((double)i)) + log(p) + logminp;
    }
    maxcp = cp[0];
    for (i = 0; i <= maxhit; i++) {
        if (maxcp < cp[i]) {
            maxcp = cp[i];
        }
    }
    for (i = 0; i <= maxhit; i++) cp[i] = exp(cp[i] - maxcp);
    for (i = 0; i <= maxhit; i++) normalize += cp[i];
    for (i = 0; i <= maxhit; i++) cp[i] /= normalize;


}

void computeCompoundPoissonDistributionKempSingleStranded(double lambda,
        int maxhit, int maxclump, double *theta, double *cp) {
    int i, j;
    double p;
    double normalize = 0.0;
    double logminp, maxcp;

#ifdef DEBUG
    Rprintf( "lambda=%f\n", lambda);
#endif
    cp[0] = -lambda;

    // compute the compound poisson distribution according to
    // Kemp in log-space
    for (i = 1; i <= maxhit; i++) {
        p = 0.0;
        j = i - maxclump + 1;
        j = (j > 0) ? j : 0;
        // find the smallest logp
        logminp = cp[j];
        for (; j < i; j++) {
            if (logminp > cp[j]) {
                logminp = cp[j];
            }
        }
        j = i - maxclump + 1;
        j = (j > 0) ? j : 0;
        // compute p which is divided by logminp
        for (; j < i; j++) {
            p += (i - j) * theta[i - j - 1] * exp(cp[j] - logminp);
        }
        // compute log(cp)
        cp[i] = log(lambda / ((double)i)) + log(p) + logminp;
    }
    // determine maxcp which is still in log-space
    maxcp = cp[0];
    for (i = 0; i <= maxhit; i++) {
        if (maxcp < cp[i]) {
            maxcp = cp[i];
        }
    }
    for (i = 0; i <= maxhit; i++) cp[i] = exp(cp[i] - maxcp);
    for (i = 0; i <= maxhit; i++) normalize += cp[i];
    for (i = 0; i <= maxhit; i++) cp[i] /= normalize;
}

// This function computes the clump size distribution
// according to Kopp et al. for both DNA strands
void clumpsizeBeta(double *beta, double *beta3p, double *beta5p,
                        double *dist, int *maxclump, int *motiflen) {

    double *theta, extention[3];
    double *delta, *deltap;

    delta = (double*)R_alloc((size_t)motiflen[0], sizeof(double));
    deltap = (double*)R_alloc((size_t)motiflen[0], sizeof(double));

    memset(delta, 0, motiflen[0]*sizeof(double));
    memset(deltap, 0, motiflen[0]*sizeof(double));

    // initialize the extention factors
    memset(extention, 0, 3 * sizeof(double));

    beta3p[0] = (beta3p[0] + EPSILON) / (1 + 2 * EPSILON);
    computeDeltas(delta, deltap, beta, beta3p, beta5p, motiflen[0]);

    computeExtentionFactorsKopp(extention, delta, deltap, beta,
                                beta3p, beta5p, motiflen[0]);
    //theta = dist;

    computeInitialClumpKopp(dist, beta3p, delta, deltap, motiflen[0]);
    computeTheta(maxclump[0], dist, extention, motiflen[0]);

}

// This function computes the clump size distribution
// according to Kopp et al. for both DNA strands
void clumpsizeBeta_singlestranded(double *beta,
                        double *dist, int *maxclump, int *motiflen) {

    double *theta, extention;
    double *delta;

    delta = (double*)R_alloc((size_t)motiflen[0], sizeof(double));

    memset(delta, 0, motiflen[0]*sizeof(double));

    // initialize the extention factors
    extention = 0.0;

    //computeDeltas(delta, deltap, beta, beta3p, beta5p, motiflen[0]);
    computeDeltasSingleStranded(delta, beta, motiflen[0]);

    //computeExtentionFactorsKopp(extention, delta, deltap, beta,
                                //beta3p, beta5p, motiflen[0]);
    computeExtentionFactorsKoppSingleStranded(&extention, beta, motiflen[0]);
    //theta = dist;
        //theta = initThetaSingleStranded(maxclumpsize);

    //computeInitialClumpKopp(theta, beta3p, delta, deltap, motiflen[0]);
    computeInitialClumpKoppSingleStranded(dist, delta, motiflen[0]);
    //computeTheta(maxclumpsize, theta, extention, motiflen[0]);
    computeThetaSingleStranded(*maxclump, dist, &extention, motiflen[0]);

}

// This function computes the clump size distribution
// according to Pape et al. for both DNA strands
void clumpsizeGamma(double *gamma, double *dist, int *maxclump, int *motiflen) {

    double *theta, extention[3];

    // initialize the extention factors
    memset(extention, 0, 3 * sizeof(double));
    gamma[motiflen[0]] = (gamma[motiflen[0]] + EPSILON) / (1 + 2 * EPSILON);
    computeExtentionFactorsPape(extention, gamma, motiflen[0]);

    //theta = dist;

    computeInitialClump(dist, gamma, motiflen[0]);
    computeTheta(*maxclump, dist, extention, motiflen[0]);
}
