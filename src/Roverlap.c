#include <R.h>
#include "overlap.h"
#include "score2d.h"


extern double Rgran, Rsiglevel;

void Roverlap(double *pfm_, int *nrow, int *ncol,
        double *alpha, double *beta, double *beta3p, double *beta5p,
        double *gamma, double *station, double *trans, int *order) {

    int i;
    double dx, pvalue;
    DMatrix pfm, cpfm;

    if (!beta||!beta3p||!beta5p) {
        error("parameters are null");
        return;
    }
    if (Rgran==0.0 || Rsiglevel==0.0) {
        error("call mdist.option  first");
        return;
    }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    cpfm.data=Calloc(nrow[0]*ncol[0],double);
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    cpfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    cpfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));
    for (i=1; i<=nrow[0]*ncol[0];i++) {
        cpfm.data[i-1]=pfm.data[nrow[0]*ncol[0]-i];
    }


    dx=(double)Rgran;
    pvalue=(double)Rsiglevel;

    computeConditionalOverlappingProbabilities(&pfm, &cpfm,
            station, trans, NULL, &pvalue, NULL, &dx, gamma, order[0]);

    for (i=1;i<pfm.nrow; i++) {
        gamma[i]/=gamma[0];
    }
    for (i=0;i<pfm.nrow; i++) {
        gamma[pfm.nrow+i]/=gamma[0];
    }
    for (i=0;i<pfm.nrow; i++) {
        gamma[pfm.nrow*2+i]/=gamma[0];
    }

    computeBetas(beta, beta3p,beta5p,gamma,pfm.nrow, 0.0);
    *alpha=gamma[0];

    Free(pfm.data);
    Free(cpfm.data);
}

void RoverlapSingleStranded(double *pfm_, int *nrow, int *ncol,
        double *alpha, double *beta, double *beta3p,
        double *beta5p, double *gamma, double *station, double *trans,
        int *order) {

    int i;
    double dx, pvalue;
    DMatrix pfm, cpfm;

    if (!beta||!beta3p||!beta5p) {
        error("parameters are null");
        return;
    }
    if (Rgran==0.0 || Rsiglevel==0.0) {
        error("call mdist.option  first");
        return;
    }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    cpfm.data=Calloc(nrow[0]*ncol[0],double);
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    cpfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    cpfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));
    for (i=1; i<=nrow[0]*ncol[0];i++) {
        cpfm.data[i-1]=pfm.data[nrow[0]*ncol[0]-i];
    }

    dx=(double)Rgran;
    pvalue=(double)Rsiglevel;

    computeConditionalOverlappingProbabilities(&pfm, &cpfm,
            station, trans, NULL, &pvalue, NULL, &dx, gamma, order[0]);

    for (i=1;i<pfm.nrow; i++) {
        gamma[i]/=gamma[0];
    }
    for (i=0;i<pfm.nrow; i++) {
        gamma[pfm.nrow+i]/=gamma[0];
    }
    for (i=0;i<pfm.nrow; i++) {
        gamma[pfm.nrow*2+i]/=gamma[0];
    }

    computeBetasSingleStranded(beta, gamma,pfm.nrow, 0.0);
    *alpha=gamma[0];

    Free(pfm.data);
    Free(cpfm.data);
}
