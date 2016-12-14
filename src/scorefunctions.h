
#ifndef scorefunctions_h
#define scorefunctions_h

#include "matrix.h"
#include "minmaxscore.h"



typedef struct {
    int length;
    int zero;
    double dx;
    int xmax, xmin;
    double (*prob)(double,double);
    double (*probinit)(double,double*, int, int);
} ScoreMetaInfo;

int initScoreMetaInfo (int min, int max, int intervalsize, double siglen, ScoreMetaInfo *meta);
// initialize the score distribution structure

void getScoresIndex(double *P,double *Q, int *score, double *dx);
double ProbBg (double p, double );
double ProbPWM (double p, double);
double ProbinitBg (double p, double *b, int, int);
double ProbinitPWM (double p, double *b, int, int);
int getScoreIndex(double P,double Q, double dx);
void getScoresInitialIndex(double *P,double *Q, int *score, double *dx, int order);
void scoreSequence(double *station, double *trans,
  DMatrix *pwm, const char *seq, int seqlen, double *scores,
  double granularity, int order);
void scoreHistogram(double *station, double *trans,
    DMatrix *pwm, const char *seq, int seqlen,
   double *dist, double granularity, int smin, int order);
#endif
