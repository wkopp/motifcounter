#include <R.h>
#include "score1d.h"

extern int Rorder;
extern double *Rstation, *Rtrans;
//extern DMatrix *&pfm, *&cpfm;
extern double Rsiglevel, Rgran;

void Rscoredist(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob) {
    MotifScore1d null;
    ExtremalScore fescore,rescore;
    double dx;
    //double pcomp;
    int i, intervalsize,maxs,mins, fmins,fmaxs,rmins, rmaxs;
    //int threshold;
    DMatrix pfm, cpfm;


    if (!Rstation||!Rtrans) {
        error("load forground and background properly");
        return;
    }
    if (!score||!prob) {
        error("parameters are null");
        return;
    }
    if (Rgran<=1e-10) { error("set granularity first! E.g. 0.1"); }

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
    //p=Rsiglevel;

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
    mins=(fmins>rmins) ? fmins : rmins;

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);
    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,Rtrans,&null, Rorder);

    computeScoreDistribution1d(&pfm, Rtrans,  Rstation,
            &null, &fescore, Rorder);

    //quantile=getQuantileWithIndex1d(&null,
                    //getQuantileIndex1d(&null.totalScore,p));
    //threshold=(int)(quantile/dx);
    //Rprintf("threshold=%e\n",quantile);
    //pcomp=getProbWithIndex1d(&null,threshold);

    for (i=0; i<null.meta.xmax-null.meta.xmin + 1; i++) {
        score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
        prob[i]=null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);
    deleteScoreDistribution1d(&null, Rorder);

    Free(pfm.data);
    Free(cpfm.data);
}

void Rscoredist_bf(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob) {
    int i;
    MotifScore1d null;
    ExtremalScore fescore,rescore;
    double dx;
    int fmins,fmaxs,rmins,rmaxs;
    int mins,maxs;
    int intervalsize;
    DMatrix pfm, cpfm;

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
    //p=Rsiglevel;

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
    mins=(fmins>rmins) ? fmins : rmins;

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);

    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,Rtrans,&null, Rorder);

    computeMarginalScoreDistribution1dBruteForce(&pfm, Rtrans,
            Rstation, &null, null.meta.xmin, Rorder);


    ////quantile=getQuantileWithIndex1d(&null,
        //getQuantileIndex1d(&null.totalScore,p));
    //threshold=(int)(quantile/dx);
    //pcomp=getProbWithIndex1d(&null,threshold);


    for (i=0; i<null.meta.xmax-null.meta.xmin+1&& i<null.meta.length; i++) {
        score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
        prob[i]=null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);
    deleteScoreDistribution1d(&null, Rorder);

    Free(pfm.data);
    Free(cpfm.data);
}
