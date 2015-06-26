
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

int numNucleotides(FILE *f) {
	int cnt=0, stop=0;
	char dig;
	int fposition=ftell(f);
	fseek(f,ftell(f)-1,SEEK_SET);
  while((dig=fgetc(f))) {
  	if (dig==EOF || dig=='>') break;
  	switch (dig) {
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'n':
		case 'N':
      cnt++;
      Rprintf("%c",dig);
      break;
    case '\n':
			break;
		case '>':
			stop=1;
  	default:
  		
  		if (dig==EOF) {
  			Rprintf("eof should have stopped!\n");
  		} else {
  		  Rprintf("something more: %c\n",dig==EOF);
  		}
  		stop=1;
  		break;
  	}
  	if(stop==1) break;
  }
  Rprintf("\n");

  fseek(f,ftell(f)-1,SEEK_SET);
  return cnt;
}

void getSeqlen(FILE *f, int *seqlen, int *numofseq) {
  char dig;
  int writeheader=0;
  int refslen=0;
  seqlen[0]=0;numofseq[0]=0;

  if (fgetc(f)!='>') {
  	error("This is not a fasta file!");
  	return;
  }

  rewind(f);

	// get the length of the first sequence
  while((dig=fgetc(f))!=EOF) {
    if (dig=='>') {
      writeheader=1;
      continue;
    }
    if (writeheader==1 && dig=='\n') {
      writeheader=0;
      continue;
    }
    if (writeheader==0) {
      refslen=numNucleotides(f);
      break;
    }
  }
  Rprintf("reflen=%d\n",refslen);

  rewind(f);

	// compare all sequences in the multiple fasta file against the
	// length of the first sequence. They must be equally long.
  while((dig=fgetc(f))) {
  	if(dig==EOF) break;
    if (dig=='>') {
      numofseq[0]++;
      writeheader=1;
    }
    if (writeheader==1 && dig=='\n') {
      writeheader=0;
      continue;
    }
    if (writeheader==0) {
    	seqlen[0]=numNucleotides(f);
    	writeheader=1;
    	if (refslen!=seqlen[0]) {
    		if (dig==EOF) Rprintf("eof observed!: ");
    		if (dig=='\n') Rprintf("\\n observed!: ");
    		error("For multiple fasta files, the individual sequences must be of"
    				" equal length. %d!=%d:dig=%c",refslen,seqlen[0],dig);
    	}
    }
  }
}

void getNucleotideSequence(FILE *f, char **seq, int seqlen, int numofseq) {
  char dig;
  int writeseq=0, writeheader=0;
  int iseqlen=0, inos=-1;



  rewind(f);
  while((dig=fgetc(f))!=EOF && inos<numofseq) {
    if (dig=='>') {
      writeheader=1;
      writeseq=0;
      inos++;
      iseqlen=0;
    }

    if (writeseq==1 && dig=='\n')  continue;

    if (writeheader==1 && dig=='\n') {
      writeheader=0;
      writeseq=1;
      continue;
    }

    if (writeseq==1 && iseqlen<seqlen) {
      seq[inos][iseqlen++]=dig;
    }
  }
}

void RnumberOfHits(char **inputfile, int *numofhits, int *slen, int *nos) {
  int Nhits,seqlen, numofseq, intervalsize;
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
  if (!numofhits||!inputfile||!slen||!nos) {
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
  getSeqlen(f, &seqlen, &numofseq);
  slen[0]=seqlen;
  nos[0]=numofseq;
  #ifdef IN_R
  seq=Calloc(numofseq, char*);
  for (i=0; i<numofseq; i++) {
    seq[i]=Calloc(seqlen+1, char);
  }
  #else
  seq=calloc(numofseq, sizeof(char*));
  for (i=0; i<numofseq; i++) {
    seq[i]=calloc(seqlen+1, sizeof(char));
  }
  #endif
  getNucleotideSequence(f, seq, seqlen, numofseq);
  fclose(f);
  //Rprintf("numofseq=%d, seqlen=%d\n", numofseq, seqlen);

  for (s=0, Nhits=0;s<numofseq; s++) {
    Nhits+=countOccurances(Rstation, Rtrans, Rpwm, Rcpwm, seq[s], seqlen, threshold, dx, Rorder);
  }

#ifdef IN_R
  for (i=0; i<numofseq; i++) Free(seq[i]);
  Free(seq);
  #else
  for (i=0; i<numofseq; i++) free(seq[i]);
  free(seq);
  #endif
  numofhits[0]=Nhits;
}

