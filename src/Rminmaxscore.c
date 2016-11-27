#include <R.h>
#include "minmaxscore.h"

extern int Rorder;
extern double *Rstation, *Rtrans;
extern double Rgran;

void Rscorerange(double *pfm_, int *nrow, int *ncol, int *scorerange) {
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs;
    double dx;
    int i;
    DMatrix pfm, cpfm;

    if (scorerange==NULL) {
        error("scorerange is null");
    }
    if (Rstation==NULL || Rtrans==NULL) {
        error("load forground and background before!");
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
    initExtremalScore(&fescore, dx, pfm.nrow, Rorder);
    initExtremalScore(&rescore, dx, cpfm.nrow, Rorder);

    loadMinMaxScores(&pfm, Rstation, Rtrans, &fescore);
    loadMinMaxScores(&cpfm, Rstation, Rtrans, &rescore);
    loadIntervalSize(&fescore, NULL);
    loadIntervalSize(&rescore, NULL);

    fmins=getTotalScoreLowerBound(&fescore);
    rmins=getTotalScoreLowerBound(&rescore);
    fmaxs=getTotalScoreUpperBound(&fescore);
    rmaxs=getTotalScoreUpperBound(&rescore);
    maxs=(fmaxs>rmaxs) ? fmaxs : rmaxs;
    mins=(fmins<rmins) ? fmins : rmins;

    scorerange[0]=maxs-mins+1;

    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);

    Free(pfm.data);
    Free(cpfm.data);
}

