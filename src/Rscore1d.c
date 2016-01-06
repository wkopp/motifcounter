#include <R.h>
#include "score1d.h"

extern int Rorder;
extern double *Rstation, *Rtrans;
extern DMatrix *Rpwm, *Rcpwm;
extern double Rsiglevel, Rgran;

void Rscoredist( double *score, double *prob) {
  MotifScore1d null, alter;
  ExtremalScore escore;
  double p, quantile, dx, pcomp;
  int i, intervalsize;
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

  initExtremalScore(&escore, dx, Rpwm->nrow, Rorder);

  loadMinMaxScores(Rpwm, Rstation, Rtrans, &escore);
  loadIntervalSize(&escore, NULL);
  //intervalsize=maxScoreIntervalSize(&escore);
  intervalsize=getTotalScoreUpperBound(&escore)-
  		getTotalScoreLowerBound(&escore)+1;

  initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
  initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &alter.meta);
  null.meta.prob=&ProbBg;
  alter.meta.prob=&ProbPWM;
  null.meta.probinit=&ProbinitBg;
  alter.meta.probinit=&ProbinitPWM;
  initScoreDistribution1d(Rpwm,Rtrans,&null, Rorder);
  initScoreDistribution1d(Rpwm,Rtrans,&alter, Rorder);

  computeScoreDistribution1d(Rpwm, Rtrans,  Rstation, 
          &null, &escore, Rorder);

  computeScoreDistribution1d(Rpwm, Rtrans,  Rstation, 
          &alter, &escore, Rorder);

  quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
  threshold=(int)(quantile/dx);
  pcomp=getProbWithIndex1d(&null,threshold);
  Rprintf("%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
        getQuantileIndex1d(&null.totalScore,p),quantile);

  for (i=0; i<=null.meta.xmax-null.meta.xmin&& i<null.meta.length; i++) {
    score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
    prob[i]=null.totalScore.y[i];
  }
  #ifdef WK
        if (output_dist!=NULL) {
          f=fopen(output_dist,"w");
          storeScoreDist1d(f,&null, 1);
          storeScoreDist1d(f,&alter, 0);
        }

  loadIntervalSize(&escore, &threshold);

  cutScoreRangeWithThreshold(&null, &escore, Rorder);
        Rprintf( "P(hit)= %f\n", getProbability1d(&null.totalScore,
          &null.meta));

        if (output_dist!=NULL) {
          fclose(f);
        }
        #endif
  deleteExtremalScore(&escore);
  deleteScoreDistribution1d(&null, Rorder);
  deleteScoreDistribution1d(&alter, Rorder);

}

void Rscoredist_bf( double *score, double *prob) {
//void Rscoredist_bruteforce(char *bgfile, char *pwmfile, char *output_quantile, 
//   char *output_dist, char *pv, char *th, char *gr) {
  int i;
  MotifScore1d null, alter;
  ExtremalScore escore;
  double p, quantile, dx, pcomp;
  int intervalsize;
  int threshold;

   dx=Rgran;
   p=Rsiglevel;


   initExtremalScore(&escore, dx, Rpwm->nrow, Rorder);

   loadMinMaxScores(Rpwm, Rstation, Rtrans, &escore);
   loadIntervalSize(&escore, NULL);
 //  intervalsize=maxScoreIntervalSize(&escore);
  intervalsize=getTotalScoreUpperBound(&escore)-
  		getTotalScoreLowerBound(&escore)+1;

   initScoreMetaInfo(getTotalScoreLowerBound(&escore),
   getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
   initScoreMetaInfo(getTotalScoreLowerBound(&escore),
   getTotalScoreUpperBound(&escore),intervalsize,dx, &alter.meta);
   null.meta.prob=&ProbBg;
   alter.meta.prob=&ProbPWM;
   null.meta.probinit=&ProbinitBg;
   alter.meta.probinit=&ProbinitPWM;
   initScoreDistribution1d(Rpwm,Rtrans,&null, Rorder);
   initScoreDistribution1d(Rpwm,Rtrans,&alter, Rorder);

   computeMarginalScoreDistribution1dBruteForce(Rpwm, Rtrans,
          Rstation, &null, null.meta.xmin, Rorder);


   quantile=getQuantileWithIndex1d(&null,
      getQuantileIndex1d(&null.totalScore,p));
   threshold=(int)(quantile/dx);
        //threshold=203;
   pcomp=getProbWithIndex1d(&null,threshold);
   Rprintf("%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
   getQuantileIndex1d(&null.totalScore,p),quantile);


  for (i=0; i<=null.meta.xmax-null.meta.xmin&& i<null.meta.length; i++) {
    score[i]=(double)(null.meta.xmin+i)*null.meta.dx;
    prob[i]=null.totalScore.y[i];
  }
  deleteExtremalScore(&escore);
  deleteScoreDistribution1d(&null, Rorder);
  deleteScoreDistribution1d(&alter, Rorder);

}

