#include <R.h>
#include "minmaxscore.h"


extern double Rgran;

void Rscorerange(double *pfm_, int *nrow, int *ncol, int *scorerange,
    double *station, double *trans, int *order) {
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs;
    double dx;
    int i;
    DMatrix pfm, cpfm;

    if (scorerange==NULL) {
        error("scorerange is null");
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
    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);
    initExtremalScore(&rescore, dx, cpfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadMinMaxScores(&cpfm, station, trans, &rescore);
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
