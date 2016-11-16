#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "forground.h"
#include "background.h"
#include "countdist.h"
#include "sequence.h"
#include "scorefunctions.h"
#include "minmaxscore.h"
#include "score1d.h"

char sampleNucleotide(double *prob) {
    double ip=unif_rand();
    double p;
    double cumprob=0.0;
    char n;
    int i;

    p=(double)(ip);
    for (i=0; i<ALPHABETSIZE; i++) {
        cumprob+=prob[i];
        if(p<=cumprob) {
            break;
        }
    }
    if (i==4)
        n= getNuc(i-1);
    else
        n= getNuc(i);
    return n;
}

void sampleInitialNucleotide(double *prob, char *seq, int order) {
    double ip=unif_rand();
    double p;
    double cumprob=0.0;
    int i, ass[order];

    p=(double)(ip);
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
        cumprob+=prob[i];
        if(p<=cumprob) {
            break;
        }
    }
    getAssignmentFromIndex(i, order, ass);
    for(i=0;i<order; i++) {
        seq[i]=getNuc(ass[i]);
    }
}

void generateRandomSequence(double *station, double *trans, char *seq, 
        int seqlen, int order) {
    int i, ind;

    sampleInitialNucleotide(station, seq, order);
    for (i=order; i<seqlen; i++) {
        ind=getIndexFromAssignment(&seq[i-order], order);
        seq[i]=sampleNucleotide(&trans[ind*ALPHABETSIZE]);
    }
}

void randomStatistics(char *seq, int seqlen, double *stat, double *trans) {
  int i;

  stat[getNucIndex(seq[0])]+=1.;
  for (i=1; i<seqlen; i++) {
    stat[getNucIndex(seq[i])]++;
    trans[getNucIndex(seq[i-1])*ALPHABETSIZE+getNucIndex(seq[i])]++;
  }
}
void normalizeStatistics(double *stat, double *trans) {
  int i, j, sum=0;
  for (i=0; i<power(ALPHABETSIZE, 1);i++) {
  sum=0;
    for (j=0; j<power(ALPHABETSIZE, 1); j++) {
      sum+=trans[i*ALPHABETSIZE +j];
    }
    for (j=0; j<power(ALPHABETSIZE, 1); j++) {
      trans[i*ALPHABETSIZE +j]/=sum;
    }
  }
  sum=0;
  for (i=0; i<power(ALPHABETSIZE, 1);i++) sum+=stat[i];
  for (i=0; i<power(ALPHABETSIZE, 1);i++) stat[i]/=sum;
}

int countOccurances(double *station, double *trans, 
  DMatrix *pwm, char *seq, int seqlen,
 int qalpha, double granularity, int order) {
 int count=0;
 int index;
 int s, i,j;
 int score[power(ALPHABETSIZE, order+1)];

  // if the sequence contains any N's, do not process the scores
  for (i=0; i<seqlen;i++) {
    if (getNucIndex(seq[i])<0) {
        return 0;
    }
  }

 for (i=0; i< seqlen-pwm->nrow+1; i++) {

   s=0;
   index=0;
   if (order>0) {
     getScoresInitialIndex(pwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s+=score[index];
   }
   for (j=order; j<pwm->nrow; j++) {

     index*=ALPHABETSIZE;
     index+=getNucIndex(seq[i+j]);

     s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index], granularity);

     index-=(index/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
   }
   if (s>=qalpha) {
     count++;
   }
 }
 return count;
}

void scoreOccurances(double *station, double *trans, 
  DMatrix *pwm, char *seq, int seqlen,
 double *dist, double granularity, int smin, int order) {
  int i, j;
  int s, index;
  int score[power(ALPHABETSIZE, order+1)];

  // if the sequence contains any N's, do not process the scores
  for (i=0; i<seqlen;i++) {
    if (getNucIndex(seq[i])<0) {
        return;
    }
  }

 for (i=0; i< seqlen-pwm->nrow+1; i++) {
   index=0;

   if (order>0) {
     getScoresInitialIndex(pwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s=score[index];
   } else {
     s=0;
   }

   for (j=order; j<pwm->nrow; j++) {

     index=index*ALPHABETSIZE + getNucIndex(seq[i+j]);

     s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index],granularity);

     index-=(index/power(ALPHABETSIZE,order))*power(ALPHABETSIZE,order);
   }
   dist[s-smin]++;
 }
}


