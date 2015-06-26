

#ifndef score2d_h
#define score2d_h

//#include "sharedtypes.h"
//#include "align.h"
#include "matrix.h"
#include "minmaxscore.h"
#include "scorefunctions.h"
#include "score1d.h"


typedef struct {
  double *y;
  int start1,end1, start2,end2;
} Score2d;

typedef struct {
  int mlen;
  ScoreMetaInfo meta;
  Score2d *ScoreBuffer1;
  Score2d *tmpScore;
  double cprob;
} MotifScore2d;



int initScoreDistribution (DMatrix *theta, DMatrix *bg1, MotifScore2d *result, int order);

int deleteMotifScore(MotifScore2d *m, int order);

//int deleteScore (Score *s);


// Compute the score distribution for the motif


int computeScoreDistribution2D(DMatrix *pwm1, DMatrix *pwm2,
  double *trans, double *station, MotifScore2d *mscore,
  ExtremalScore *escore1, ExtremalScore *escore2,
  int shift, int nsubmotif,
  MotifScore1d *init1d, MotifScore2d *init2d, int order);
void computeConditionalOverlappingProbabilities(DMatrix *pwm1, DMatrix *pwm2, double *station,
  double *trans, FILE *fout, double *pvalue, int *inth, double *dx, double *gamma, int order);
#ifdef WKOPP
int computeScoreDistribution2D(DMatrix *pwm1, DMatrix *pwm2,
  double *trans, double *station, MotifScore2d *mscore,
  ExtremalScore *fext, ExtremalScore *rext,
  MotifScore1d *init1d, MotifScore2d *init2d, int, int, int);
#endif
//int computeScoreDistribution2D(DMatrix *pwm, DMatrix *cpwm,
//  DMatrix *ftrans, DMatrix *rtrans, DMatrix *station, MotifScore2d *mscore);
// sums over all values of S greater or equal than threshold
//int computeConditionalScoreDistribution(MotifScore *s, double threshold);

void scoredist2d(char *bgfile, char *pwmfile, char *oq,char*od, char*,char*, char*);
void scoredist2d_bfinit(char *bgfile, char *pwmfile, 
  char *output_quantile, char *output_dist, char *pv, char *gr,
    char *sth);
void scoredist2d_bruteforce(char *bgfile, char *pwmfile, char *oq,char*od, char*gr,char*sth,char*pv);
// This function is indented to compute the conditional distribution
// over the scores, given the current score exceeds the threshold Sth.
//int computeConditionalScore(DMatrix *theta, DMatrix *bg1, DMatrix *bg0, 
// DMatrix *thresholds, int index, Signal *dist);

// This function computes a matrix with maximum achievable scores at
// the current position, given the current letter.
#endif
