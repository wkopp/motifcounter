
#include <time.h>
#include <R.h>
#include "matrix.h"
#include "simulate.h"
#include "sequence.h"
#include "minmaxscore.h"
#include "score1d.h"

extern int Rorder;
extern DMatrix *Rpwm, *Rcpwm;
extern double *Rstation, *Rtrans;
extern double Rsiglevel, Rgran;

void RnumberOfHits(char **inputfile, int *numofhits, int *nseq, int *lseq,  
        int *singlestranded) {
    int intervalsize;
    ExtremalScore escore;
    MotifScore1d null;
    Sequence seq;
    double dx, quantile,pvalue;
    int s;
    int threshold;
    FILE *f;

    if (!Rpwm||!Rcpwm||!Rstation||!Rtrans) {
        error("load forground and background properly");
        return;
    }
    if (!numofhits||!inputfile||!nseq||!lseq) {
        error("parameters are null");
        return;
    }
    if (Rsiglevel==0.0 || Rgran==0.0) {
        error ("call mdist.option first");
        return;
    }

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

    quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,
                pvalue));
    threshold=(int)(quantile/dx);
    deleteExtremalScore(&escore);
    deleteScoreDistribution1d(&null, Rorder);

    f=fopen(inputfile[0], "r");
    if (!f) {
        error("no such file: %s\n", inputfile[0]);
        return;
    }
    allocSequence(&seq,*nseq, lseq);
    getSequence(f,&seq);
    fclose(f);

    for (s=0;s<*nseq; s++) {
        numofhits[s]+=countOccurances(Rstation, Rtrans, Rpwm, seq.seq[s], 
                seq.lseq[s], threshold, dx, Rorder);
        if (*singlestranded==0) {
            numofhits[s]+=countOccurances(Rstation, Rtrans, Rcpwm, seq.seq[s], 
                    seq.lseq[s], threshold, dx, Rorder);
        }
    }

    destroySequence(&seq);
}

