#include <R.h>
#include "minmaxscore.h"
//#include "

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
    return;
  }
  if (Rpwm==NULL || Rcpwm==NULL || Rstation==NULL || Rtrans==NULL) {
    error("load forground and background before!");
    return;
  }

  dx=Rgran;
//  Rprintf("range=%d, gran=%e\n",scorerange[0],dx);
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
//	Rprintf("Rscorerange: f[lb, ub]=[%d,%d]\n",fmins, fmaxs);
//	Rprintf("Rscorerange: r[lb, ub]=[%d,%d]\n",rmins, rmaxs);
  maxs=(fmaxs>rmaxs) ? fmaxs : rmaxs;
  mins=(fmins<rmins) ? fmins : rmins;
//	Rprintf("Rscorerange: [lb, ub]=[%d,%d]\n",mins, maxs);

  scorerange[0]=maxs-mins+1;

  deleteExtremalScore(&fescore);
  deleteExtremalScore(&rescore);
 // Rprintf("range=%d\n",scorerange[0]);
}

