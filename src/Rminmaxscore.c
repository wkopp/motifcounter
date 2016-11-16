#include <R.h>
#include "minmaxscore.h"

extern DMatrix *Rpwm, *Rcpwm;
extern int Rorder;
extern double *Rstation, *Rtrans;
extern double Rgran;

void Rscorerange(int *scorerange) {
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs;
    double dx;

    if (scorerange==NULL) {
        error("scorerange is null");
    }
    if (Rpwm==NULL || Rcpwm==NULL || Rstation==NULL || Rtrans==NULL) {
        error("load forground and background before!");
    }

    dx=Rgran;
    initExtremalScore(&fescore, dx, Rpwm->nrow, Rorder);
    initExtremalScore(&rescore, dx, Rcpwm->nrow, Rorder);

    loadMinMaxScores(Rpwm, Rstation, Rtrans, &fescore);
    loadMinMaxScores(Rcpwm, Rstation, Rtrans, &rescore);
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
}

