
#ifndef scorefunctions_h
#define scorefunctions_h

//#include "sharedtypes.h"
//#include "align.h"
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

//void getScores(double *P,double *Q, double *score, double *dx);
void getScoresIndex(double *P,double *Q, int *score, double *dx);
double ProbBg (double p, double );
double ProbPWM (double p, double);
double ProbinitBg (double p, double *b, int, int);
double ProbinitPWM (double p, double *b, int, int);
double getDiscretizedScore(double P, double Q, double dx);
int getScoreIndex(double P,double Q, double dx);
void getScoresInitialIndex(double *P,double *Q, int *score, double *dx, int order);
void printSeq(int index, int len);
//int deleteScore (Score *s);



// This function is indented to compute the conditional distribution
// over the scores, given the current score exceeds the threshold Sth.
//int computeConditionalScore(DMatrix *theta, DMatrix *bg1, DMatrix *bg0, 
// DMatrix *thresholds, int index, Signal *dist);

// This function computes a matrix with maximum achievable scores at
// the current position, given the current letter.
#endif
