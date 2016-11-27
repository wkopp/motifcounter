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
extern double *Rstation, *Rtrans;
extern double *RstationForSampling, *RtransForSampling;
extern double Rsiglevel, Rgran;

void RsimulateCountDistribution(double *pfm_, int *nrow, int *ncol,
    double *distribution, int* perm, 
    int *nseq, int *lseq, int *mxhit, int *singlestranded) {
    int seqlen,maxhits,Nperm, intervalsize;
    int *_nh;
    ExtremalScore escore;
    MotifScore1d null;
    double dx, quantile,pvalue;
    int n, s;
    int i;
    char *seq;
    DMatrix pfm, cpfm;

    int threshold;

    if (!Rstation||!RstationForSampling) {
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
    GetRNGstate();

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

    maxhits=mxhit[0];
    Nperm=perm[0];
    pvalue=Rsiglevel;
    dx=Rgran;

    // compute significance threshold
    initExtremalScore(&escore, dx, pfm.nrow, Rorder);
    loadMinMaxScores(&pfm, Rstation, Rtrans, &escore);
    loadIntervalSize(&escore, NULL);

    intervalsize=getTotalScoreUpperBound(&escore)-
                    getTotalScoreLowerBound(&escore)+1;

    initScoreMetaInfo(getTotalScoreLowerBound(&escore),
            getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
    null.meta.prob=&ProbBg;
    null.meta.probinit=&ProbinitBg;
    initScoreDistribution1d(&pfm,Rtrans,&null, Rorder);

    computeScoreDistribution1d(&pfm, Rtrans,  Rstation, &null, &escore, Rorder);

    quantile=getQuantileWithIndex1d(&null,
                getQuantileIndex1d(&null.totalScore,pvalue));
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

    if (*singlestranded==0) {
        Rprintf("count in both strands\n");
    }else{
        Rprintf("count single  strands\n");
    }

    for (n=0; n<Nperm; n++) {
        _nh[n]=0;
        for (s=0;s<*nseq; s++) {
            generateRandomSequence(RstationForSampling, RtransForSampling, seq, 
                    lseq[s], RorderForSampling);
            _nh[n]+=countOccurances(Rstation, Rtrans, &pfm, seq, lseq[s], 
                    threshold, dx, Rorder);
            if (*singlestranded==0) {
                _nh[n]+=countOccurances(Rstation, Rtrans, &cpfm, seq, lseq[s], 
                        threshold, dx, Rorder);
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
    Free(pfm.data);
    Free(cpfm.data);
}

void RsimulateScores(double *pfm_, int *nrow, int *ncol,
    double *scores, double *distribution, int *slen, 
    int *perm) {
    int i, n;
    char *seq=NULL;
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs;
    int seqlen, Nperm;
    double dx;
    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }
    if (RstationForSampling==NULL || RtransForSampling==NULL) {
        error("Background model uninitialized! "
                "Use readBackgroundForSampling()");
        return;
    }

    GetRNGstate();

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

    seqlen=slen[0];
    Nperm=perm[0];

    dx=Rgran;

    seq=Calloc(seqlen+1, char);
    if(seq==NULL) {
        error("Allocation for sequence failed");
    }

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

    for (n=0; n<Nperm; n++) {
        generateRandomSequence(RstationForSampling, 
                RtransForSampling, seq, seqlen, RorderForSampling);
        scoreOccurances(Rstation, Rtrans, 
                &pfm, seq, seqlen, distribution, dx, mins,Rorder);
    }
    for (n=0; n<maxs-mins+1; n++) {
        distribution[n]/=(double)Nperm*(seqlen-pfm.nrow+1);
    }

    for (i=0; i<maxs-mins+1; i++) {
        scores[i]= (double)(mins+i)*dx;
    }
    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);
    PutRNGstate();
    Free(seq);
    Free(pfm.data);
    Free(cpfm.data);
}

