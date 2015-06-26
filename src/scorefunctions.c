#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "scorefunctions.h"
#include "sequence.h"
#include "forground.h"
//#include "inputoutput.h"

double ProbinitPWM (double b, double *f, int index, int order) {
  double g=1;
  int o, l;

  for (o=0; o<order; o++) {
    l=index/power(ALPHABETSIZE, order-o-1);
    g*=f[o*ALPHABETSIZE+l];
    index-=l*power(ALPHABETSIZE, order-o-1);
  }
  return g;
}
double ProbinitBg (double b, double *f, int ass, int order) {
  return b;
}

double ProbPWM (double b, double f) {

  return f;
}
double ProbBg (double b, double f) {
  return b;
}

void getScoresIndex(double *P,double *Q, int *score, double *dx) {
  int i;

  for (i=0;i<ALPHABETSIZE; i++) {
    score[i]=getScoreIndex(P[i],Q[i],*dx);
  }
}

void getScoresInitialIndex(double *P,double *Q, int *score, double *dx, int order) {
  int i, j;
  int ass[order];
  double s;

  if (order==0) {
    order++;
  }
  for (i=0;i<power(ALPHABETSIZE, order); i++) {
    s=0;
    getAssignmentFromIndex(i, order, ass);
    for (j=0; j<order; j++) {
      s+=log(P[j*ALPHABETSIZE+ass[j]]);
    }
    s-=log(Q[i]);
    score[i]=(int)roundl(s/ (*dx));
    //fprintf(stdout, "si=%d, s=%f\n",score[i], s);
  }
}

double getScore(double P,double Q) {
  return log(P/Q);
}

int getScoreIndex(double P,double Q, double dx) {

  return (int)roundl(getScore(P,Q)/dx);
}

double getDiscretizedScore(double P, double Q, double dx) {
  double s=(double)getScoreIndex(P,Q,dx);
  return dx* s;
}

int initScoreMetaInfo (int smin, int smax, int intervalsize, double dx, ScoreMetaInfo *meta) {
  meta->length=(intervalsize)+1;
  meta->dx=dx;
  meta->xmax=smax;
  meta->xmin=smin;
  meta->zero=0;

  meta->prob=&ProbBg;
  meta->probinit=&ProbinitBg;
  return 0;
}

void printSeq(int index, int len) {
#ifndef IN_R
  int tmp, i;
  char nuc[]="acgt";
  tmp=index;
  for (i=0; i<len; i++) {
    fprintf(stdout,"%c", nuc[tmp/power(ALPHABETSIZE,len-i-1)]);
    tmp-=(tmp/power(ALPHABETSIZE,len-i-1))*power(ALPHABETSIZE,len-i-1);
  }
  #endif
}


