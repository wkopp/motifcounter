#include <R.h>
#include "score1d.h"


extern double Rsiglevel, Rgran;

void Rscoredist(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob,double *station, double *trans, int *order) {
    MotifScore1d null;
    ExtremalScore fescore,rescore;
    double dx;
    //double pcomp;
    int i, intervalsize,maxs,mins, fmins,fmaxs,rmins, rmaxs;
    //int threshold;
    DMatrix pfm, cpfm;


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
    mins=(fmins>rmins) ? fmins : rmins;

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);
    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,trans,&null, order[0]);

    computeScoreDistribution1d(&pfm, trans,  station,
            &null, &fescore, order[0]);

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
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
    Free(cpfm.data);
}

void Rscoredist_bf(double *pfm_, int *nrow, int *ncol,
    double *score, double *prob, double *station, double *trans, int *order) {
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
    mins=(fmins>rmins) ? fmins : rmins;

    intervalsize=maxs-mins;

    initScoreMetaInfo(mins,
            maxs,intervalsize,dx, &null.meta);

    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,trans,&null, order[0]);

    computeMarginalScoreDistribution1dBruteForce(&pfm, trans,
            station, &null, null.meta.xmin, order[0]);


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
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
    Free(cpfm.data);
}
