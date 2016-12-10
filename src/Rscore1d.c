#include <R.h>
#include "score1d.h"


extern double Rsiglevel, Rgran;

void Rscoredist(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob,double *station, double *trans, int *order) {
    MotifScore1d null;
    ExtremalScore fescore;
    double dx;
    //double pcomp;
    int i, intervalsize,maxs,mins;
    //int threshold;
    DMatrix pfm;


    if (!score||!prob) {
        error("parameters are null");
        return;
    }
    if (Rgran<=1e-10) { error("set granularity first! E.g. 0.1"); }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL) {
        error("Rscoredist: Memory allocation failed");
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

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);
    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,trans,&null, order[0]);

    computeScoreDistribution1d(&pfm, trans,  station,
            &null, &fescore, order[0]);

    for (i=0; i<null.meta.xmax-null.meta.xmin + 1; i++) {
        score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
        prob[i]=null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
}

void Rscoredist_bf(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob, double *station, double *trans, int *order) {
    int i;
    MotifScore1d null;
    ExtremalScore fescore;
    double dx;
    int mins,maxs;
    int intervalsize;
    DMatrix pfm;

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL) {
        error("Rscoredist_bf: Memory allocation failed");
    }
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));

    dx=Rgran;
    //p=Rsiglevel;

    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins=getTotalScoreLowerBound(&fescore);
    maxs=getTotalScoreUpperBound(&fescore);

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);

    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,trans,&null, order[0]);

    computeMarginalScoreDistribution1dBruteForce(&pfm, trans,
            station, &null, null.meta.xmin, order[0]);

    for (i=0; i<null.meta.xmax-null.meta.xmin+1&& i<null.meta.length; i++) {
        score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
        prob[i]=null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
}
