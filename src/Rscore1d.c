#include <R.h>
#include "score1d.h"

extern int Rorder;
extern double *Rstation, *Rtrans;
extern DMatrix *Rpwm, *Rcpwm;
extern double Rsiglevel, Rgran;

void Rscoredist( double *score, double *prob) {
  MotifScore1d null;
  ExtremalScore fescore,rescore;
  double p, quantile, dx, pcomp;
  int i, intervalsize,maxs,mins, fmins,fmaxs,rmins, rmaxs;
  int threshold;


  if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
    error("load forground and background properly");
    return;
  }
  if (!score||!prob) {
    error("parameters are null");
    return;
  }
  if (Rgran<=1e-10) { error("set granularity first! E.g. 0.1"); }
  dx=Rgran;
  p=Rsiglevel;

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
  mins=(fmins>rmins) ? fmins : rmins;

  //initExtremalScore(&escore, dx, Rpwm->nrow, Rorder);

  //loadMinMaxScores(Rpwm, Rstation, Rtrans, &escore);
  //loadIntervalSize(&escore, NULL);
  //intervalsize=getTotalScoreUpperBound(&escore)-
 // 		getTotalScoreLowerBound(&escore)+1;
  intervalsize=maxs-mins;

  initScoreMetaInfo(mins,
           maxs,intervalsize,dx, &null.meta);
  null.meta.prob=&ProbBg;
  null.meta.probinit=&ProbinitBg;
  initScoreDistribution1d(Rpwm,Rtrans,&null, Rorder);

  computeScoreDistribution1d(Rpwm, Rtrans,  Rstation, 
          &null, &fescore, Rorder);

  quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
  threshold=(int)(quantile/dx);
  pcomp=getProbWithIndex1d(&null,threshold);
  Rprintf("%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
        getQuantileIndex1d(&null.totalScore,p),quantile);

  //for (i=0; i<null.meta.xmax-null.meta.xmin + 1 && i<null.meta.length; i++) {
  Rprintf("min:%d, max:%d\n",null.meta.xmin,null.meta.xmax);
  for (i=0; i<null.meta.xmax-null.meta.xmin + 1; i++) {
    score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
    prob[i]=null.totalScore.y[i];
  }
  deleteExtremalScore(&fescore);
  deleteExtremalScore(&rescore);
  deleteScoreDistribution1d(&null, Rorder);

}

void Rscoredist_bf( double *score, double *prob) {
  int i;
  MotifScore1d null;
  ExtremalScore fescore,rescore;
  double p, quantile, dx, pcomp;
  int fmins,fmaxs,rmins,rmaxs;
  int mins,maxs;
  int intervalsize;
  int threshold;

   dx=Rgran;
   p=Rsiglevel;

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
  mins=(fmins>rmins) ? fmins : rmins;

  //initExtremalScore(&escore, dx, Rpwm->nrow, Rorder);

  //loadMinMaxScores(Rpwm, Rstation, Rtrans, &escore);
  //loadIntervalSize(&escore, NULL);
  //intervalsize=getTotalScoreUpperBound(&escore)-
 // 		getTotalScoreLowerBound(&escore)+1;
  intervalsize=maxs-mins;

  initScoreMetaInfo(mins,
           maxs,intervalsize,dx, &null.meta);

   null.meta.prob=&ProbBg;
   null.meta.probinit=&ProbinitBg;
   initScoreDistribution1d(Rpwm,Rtrans,&null, Rorder);

   computeMarginalScoreDistribution1dBruteForce(Rpwm, Rtrans,
          Rstation, &null, null.meta.xmin, Rorder);


   quantile=getQuantileWithIndex1d(&null,
      getQuantileIndex1d(&null.totalScore,p));
   threshold=(int)(quantile/dx);
   pcomp=getProbWithIndex1d(&null,threshold);
   Rprintf("%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
   getQuantileIndex1d(&null.totalScore,p),quantile);


  for (i=0; i<null.meta.xmax-null.meta.xmin+1&& i<null.meta.length; i++) {
    score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
    prob[i]=null.totalScore.y[i];
  }
  deleteExtremalScore(&fescore);
  deleteExtremalScore(&rescore);
  deleteScoreDistribution1d(&null, Rorder);

}

