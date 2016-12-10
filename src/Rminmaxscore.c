#include <R.h>
#include "minmaxscore.h"


extern double Rgran;

void Rscorerange(double *pfm_, int *nrow, int *ncol, int *scorerange,
    double *station, double *trans, int *order) {
    ExtremalScore fescore;
    int mins, maxs;
    double dx;
    int i;
    DMatrix pfm;

    if (scorerange==NULL) {
        error("scorerange is null");
    }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL) {
        error("Rscorerange: memory allocation failed");
    }
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));

    dx=Rgran;
    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins=getTotalScoreLowerBound(&fescore);
    maxs=getTotalScoreUpperBound(&fescore);

    scorerange[0]=maxs-mins+1;

    deleteExtremalScore(&fescore);

    Free(pfm.data);
}
