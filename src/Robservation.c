
#include <time.h>
#include <R.h>
#include "matrix.h"
#include "simulate.h"
#include "minmaxscore.h"
#include "score1d.h"

extern int Rorder;
extern DMatrix *Rpwm, *Rcpwm;
extern double *Rstation, *Rtrans;
extern double Rsiglevel, Rgran;

void getNucleotideSequence(FILE *f, char **seq, int *nseq, int *lseq) {
  char *buffer=NULL;
  int writeseq=0, writeheader=0, i=0;
  int iseq;
  int ipos;
  int lmax=-1;

	//Rprintf("nseq=%d,lseq=%d",nseq[0],lseq[0]);
  for (i=0; i< *nseq; i++) {
    if (lmax<lseq[i]) {
      lmax=lseq[i];
    }
  }
	//Rprintf("lmax=%d\n",lmax);
  buffer=Calloc(lmax+1, char);

  iseq=0;
  ipos=0;
  while(fgets(buffer, (lmax+1)*sizeof(char), f)!=NULL) {

    //Rprintf(".%s.",buffer);
    for (i=0; i<strlen(buffer); i++) {
      if (buffer[i]=='>') {
        iseq++;
        ipos=0;
        writeheader=1;
        writeseq=0;
      }
      if (writeheader==1 && buffer[i]=='\n') { 
      	writeheader=0; writeseq=1; break; 
      }
      if (writeseq==1 && buffer[i]=='\n') break;
      if (writeseq==1 && isNucleotide(buffer[i])==1) {
      	//Rprintf("iseq=%d, ipos=%d, '%c'",iseq, ipos, buffer[i]);
      	seq[iseq-1][ipos++]= buffer[i];
      	//Rprintf("-->%c",buffer[i]);
			}
      if (writeseq==1 && isNucleotide(buffer[i])<0) {
      	seq[iseq-1][0]= 0;
        warning("Sequence number %d contains 'n' or 'N' and is discarded.",iseq);
        writeseq=0;
        break;
      }
    }
  }

  if (buffer) Free(buffer);
}

void RnumberOfHits(char **inputfile, int *numofhits, int *nseq, int *lseq) {
  int Nhits,  intervalsize;
  ExtremalScore escore;
  MotifScore1d null;
  double dx, quantile,pvalue;
  int s, i;
  char **seq;
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
  intervalsize=maxScoreIntervalSize(&escore);

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

  f=fopen(inputfile[0], "r");
  if (!f) {
    error("no such file: %s\n", inputfile[0]);
    return;
  }
  //getSeqlen(f, &seqlen, &numofseq);
  //slen[0]=seqlen;
  //nos[0]=numofseq;
  seq=Calloc(*nseq, char*);
  for (i=0; i<*nseq; i++) {
    seq[i]=Calloc(lseq[i]+1, char);
  }
  //Rprintf("until here alright!\n");

  getNucleotideSequence(f, seq, nseq, lseq);
  fclose(f);
  //Rprintf("until here alright!\n");
  //Rprintf("numofseq=%d, seqlen=%d\n", numofseq, seqlen);

  for (s=0, Nhits=0;s<*nseq; s++) {
    Nhits+=countOccurances(Rstation, Rtrans, Rpwm, Rcpwm, seq[s], lseq[s], threshold, dx, Rorder);
  }

  for (i=0; i<*nseq; i++) Free(seq[i]);
  Free(seq);
  numofhits[0]=Nhits;
}

void RnumberOfHitsSingleStranded(char **inputfile, int *numofhits, int *nseq, int *lseq) {
  int Nhits,  intervalsize;
  ExtremalScore escore;
  MotifScore1d null;
  double dx, quantile,pvalue;
  int s, i;
  char **seq;
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
  intervalsize=maxScoreIntervalSize(&escore);

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

  f=fopen(inputfile[0], "r");
  if (!f) {
    error("no such file: %s\n", inputfile[0]);
    return;
  }
  seq=Calloc(*nseq, char*);
  for (i=0; i<*nseq; i++) {
    seq[i]=Calloc(lseq[i]+1, char);
  }

  getNucleotideSequence(f, seq, nseq, lseq);
  fclose(f);
  //Rprintf("numofseq=%d, seqlen=%d\n", numofseq, seqlen);

  for (s=0, Nhits=0;s<*nseq; s++) {
    Nhits+=countOccurancesSingleStranded(Rstation, 
    		Rtrans, Rpwm, Rcpwm, seq[s], lseq[s], threshold, dx, Rorder);
  }

  for (i=0; i<*nseq; i++) Free(seq[i]);
  Free(seq);
  numofhits[0]=Nhits;
}

