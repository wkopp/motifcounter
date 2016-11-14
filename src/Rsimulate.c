#include <time.h>
#include <R.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h"
#include "simulate.h"
#include "minmaxscore.h"
#include "background.h"
#include "score1d.h"

extern int Rorder;
extern int RorderForSampling;
extern DMatrix *Rpwm, *Rcpwm;
extern double *Rstation, *Rtrans;
extern double *RstationForSampling, *RtransForSampling;
extern double Rsiglevel, Rgran;

void RsimulateCountDistribution( double *distribution, int* perm, 
   int *nseq, int *lseq, int *mxhit, int *singlestranded) {
  int seqlen,maxhits,Nperm, intervalsize;
  int *_nh;
  ExtremalScore escore;
  MotifScore1d null;
  double dx, quantile,pvalue;
  int n, s;
  int i;
  char *seq;

  int threshold;

  if (!Rpwm||!Rcpwm||!Rstation||!RstationForSampling) {
    error("load forground and background properly");
    return;
  }
  if (!distribution||!perm||!nseq||!lseq||!mxhit) {
    error("parameters are null");
    return;
  }
  if (Rgran==0.0 || Rsiglevel==0.0) {
    error("call mdist.option  first");
    return;
  }
  Rprintf("perm=%d, nseq=%d, lseq=%d,maxhits=%d,singlestranded=%d\n",
          perm[0],nseq[0],lseq[0],mxhit[0],singlestranded[0]);
  GetRNGstate();

  maxhits=mxhit[0];
  Nperm=perm[0];
  pvalue=Rsiglevel;
  dx=Rgran;

  // compute significance threshold
  initExtremalScore(&escore, dx, Rpwm->nrow, Rorder);
  loadMinMaxScores(Rpwm, Rstation, Rtrans, &escore);
  loadIntervalSize(&escore, NULL);

  intervalsize=getTotalScoreUpperBound(&escore)-
		getTotalScoreLowerBound(&escore)+1;

  initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
  null.meta.prob=&ProbBg;
  null.meta.probinit=&ProbinitBg;
  initScoreDistribution1d(Rpwm,Rtrans,&null, Rorder);

  computeScoreDistribution1d(Rpwm, Rtrans,  Rstation, &null, &escore, Rorder);

  quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,pvalue));
  threshold=(int)(quantile/dx);
  deleteExtremalScore(&escore);
  deleteScoreDistribution1d(&null, Rorder);
  // end  of significance threshold computation

  seqlen=0;
  for (i=0; i<*nseq; i++) {
    if(seqlen<lseq[i]) {
      seqlen=lseq[i];
    }
  }
  seq=Calloc(seqlen+1,char);
  if (seq==NULL) {
  	error("Memory-allocation in RsimulateCountDistribution failed");
	}
  _nh=Calloc(Nperm,int);
  if (_nh==NULL) {
  	error("Memory-allocation in RsimulateCountDistribution failed");
	}

  Rprintf("singlestranded=%d\n",*singlestranded);
  if (*singlestranded==0) {
      Rprintf("count in both strands\n");
  }else{
      Rprintf("count single  strands\n");
  }

  for (n=0; n<Nperm; n++) {
    _nh[n]=0;
    for (s=0;s<*nseq; s++) {
      generateRandomSequence(RstationForSampling, RtransForSampling, seq, lseq[s], RorderForSampling);
      _nh[n]+=countOccurances(Rstation, Rtrans, Rpwm, seq, lseq[s], threshold, dx, Rorder);
      if (*singlestranded==0) {
        _nh[n]+=countOccurances(Rstation, Rtrans, Rcpwm, seq, lseq[s], threshold, dx, Rorder);
      }

    }
    if (_nh[n]>maxhits) _nh[n]=maxhits;
  }
  for (n=0; n<Nperm; n++) {
    distribution[_nh[n]]+=1.0/(double)Nperm;
  }

  Free(_nh);
  Free(seq);
  PutRNGstate();
}

void RsimulateScores(double *scores, double *distribution, int *slen, 
  int *perm) {
  int i, n;
  char *seq=NULL;
  ExtremalScore fescore, rescore;
  int fmins,fmaxs, rmins,rmaxs;
  int mins, maxs;
  int seqlen, Nperm;
  double dx;
  //double *rstat, *rtrans;

  if (Rgran==0.0) {
    error("call mdistOption  first");
    return;
  }
  if (RstationForSampling==NULL || RtransForSampling==NULL) {
    error("Background model uninitialized! Use readBackgroundForSampling()");
    return;
  }

  //srand(time(0));
  GetRNGstate();

  seqlen=slen[0];
  Nperm=perm[0];

  dx=Rgran;
  Rprintf("len=%d, perm=%d, gran=%e\n", seqlen, Nperm, dx);

  seq=Calloc(seqlen+1, char);
  if(seq==NULL) {
   	error("Allocation for sequence failed");
  }

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

  for (n=0; n<Nperm; n++) {
    generateRandomSequence(RstationForSampling, 
            RtransForSampling, seq, seqlen, RorderForSampling);
    scoreOccurances(Rstation, Rtrans, 
            Rpwm, seq, seqlen, distribution, dx, mins,Rorder);
  }
  for (n=0; n<maxs-mins+1; n++) {
    distribution[n]/=(double)Nperm*(seqlen-Rpwm->nrow+1);
  }

  Rprintf("min:%d, max:%d\n",mins,maxs);
  for (i=0; i<maxs-mins+1; i++) {
    scores[i]= (double)(mins+i)*dx;
  }
  deleteExtremalScore(&fescore);
  deleteExtremalScore(&rescore);
  PutRNGstate();
  Free(seq);
}

