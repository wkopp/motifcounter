#include <R.h>
#include "sequence.h"
#include "matrix.h"
#include "scorefunctions.h"

extern double Rgran;

void Rscoresequence(double *pfm_, int *nrow, int *ncol, char **seq,
    double *fscores, double *rscores, int *slen,
    double *station, double *trans, int *order) {
    int i;

    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }


    pfm.data=Calloc(nrow[0]*ncol[0],double);
    cpfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL || cpfm.data==NULL) {
        error("Rscoresequence: Memory allocation failed");
    }
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    cpfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    cpfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));
    for (i=1; i<=nrow[0]*ncol[0];i++) {
        cpfm.data[i-1]=pfm.data[nrow[0]*ncol[0]-i];
    }

    scoreSequence(station, trans,
        &pfm, seq[0], slen[0], fscores,
        Rgran, order[0]);
    scoreSequence(station, trans,
        &cpfm, seq[0], slen[0], rscores,
        Rgran, order[0]);

    Free(pfm.data);
    Free(cpfm.data);
}


void RscoreHistogram(double *pfm_, int *nrow, int *ncol, char **seq, int *slen,
    double *scorebins,  double *frequency,
    double *station, double *trans, int *order) {
    int i;
    ExtremalScore fescore;
    int mins, maxs, noscores;
    DMatrix pfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }


    pfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL) {
        error("RscoreHistogram: Memory allocation failed");
    }
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));

    initExtremalScore(&fescore, Rgran, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins=getTotalScoreLowerBound(&fescore);
    maxs=getTotalScoreUpperBound(&fescore);

    // if the sequence contains any N's, do not process the scores
    noscores=0;
    for (i=0; i<slen[0];i++) {
        if (getNucIndex(seq[0][i])<0) {
            noscores=1;
            break;
        }
    }
    if (noscores==0) {
        scoreHistogram(station, trans,
            &pfm, seq[0], slen[0], frequency, Rgran, mins,order[0]);
    }
    for (i=0; i<maxs-mins+1; i++) {
        scorebins[i]= (double)(mins+i)*Rgran;
    }
    deleteExtremalScore(&fescore);

    Free(pfm.data);
}
