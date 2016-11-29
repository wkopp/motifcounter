#include <R.h>
#include "sequence.h"
#include "matrix.h"
#include "scorefunctions.h"

extern double Rgran;
extern double *Rstation, *Rtrans;
extern int Rorder;


void Rscoresequence(double *pfm_, int *nrow, int *ncol, char **seq,
    double *fscores, double *rscores, int *slen) {
    int i, n;

    double dx;
    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }
    if (Rstation==NULL || Rtrans==NULL) {
        error("Background model uninitialized! "
                "Use readBackground()");
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

    dx=Rgran;
    scoreSequence(Rstation, Rtrans,
        &pfm, seq[0], slen[0], fscores,
        Rgran, Rorder);
    scoreSequence(Rstation, Rtrans,
        &cpfm, seq[0], slen[0], rscores,
        Rgran, Rorder);

    Free(pfm.data);
    Free(cpfm.data);
}


void RscoreHistogram(double *pfm_, int *nrow, int *ncol, char **seq, int *slen,
    double *scorebins,  double *frequency) {
    int i, n;
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs, noscores;
    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }
    if (Rstation==NULL || Rtrans==NULL) {
        error("Background model uninitialized! "
                "Use readBackgroundForSampling()");
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

    initExtremalScore(&fescore, Rgran, pfm.nrow, Rorder);
    initExtremalScore(&rescore, Rgran, cpfm.nrow, Rorder);

    loadMinMaxScores(&pfm, Rstation, Rtrans, &fescore);
    loadMinMaxScores(&cpfm, Rstation, Rtrans, &rescore);
    loadIntervalSize(&fescore, NULL);
    loadIntervalSize(&rescore, NULL);

    fmins=getTotalScoreLowerBound(&fescore);
    rmins=getTotalScoreLowerBound(&rescore);
    fmaxs=getTotalScoreUpperBound(&fescore);
    rmaxs=getTotalScoreUpperBound(&rescore);
    maxs=(fmaxs>rmaxs) ? fmaxs : rmaxs;
    mins=(fmins>rmins) ? fmins : rmins;

    // if the sequence contains any N's, do not process the scores
    noscores=0;
    for (i=0; i<slen[0];i++) {
        if (getNucIndex(seq[0][i])<0) {
            noscores=1;
            break;
        }
    }
    if (noscores==0) {
        scoreHistogram(Rstation, Rtrans,
            &pfm, seq[0], slen[0], frequency, Rgran, mins,Rorder);
    }
    for (i=0; i<maxs-mins+1; i++) {
        scorebins[i]= (double)(mins+i)*Rgran;
    }
    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);

    Free(pfm.data);
    Free(cpfm.data);
}
