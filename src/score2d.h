

#ifndef score2d_h
#define score2d_h

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

// Compute the score distribution for the motif


int computeScoreDistribution2D(DMatrix *pwm1, DMatrix *pwm2,
  double *trans, double *station, MotifScore2d *mscore,
  ExtremalScore *escore1, ExtremalScore *escore2,
  int shift, int nsubmotif,
  MotifScore1d *init1d, MotifScore2d *init2d, int order);
void computeConditionalOverlappingProbabilities(DMatrix *pwm1, DMatrix *pwm2, double *station,
  double *trans, FILE *fout, double *pvalue, int *inth, double *dx, double *gamma, int order);

void scoredist2d(char *bgfile, char *pwmfile, char *oq,char*od, char*,char*, char*);
void scoredist2d_bfinit(char *bgfile, char *pwmfile, 
  char *output_quantile, char *output_dist, char *pv, char *gr,
    char *sth);
void scoredist2d_bruteforce(char *bgfile, char *pwmfile, char *oq,char*od, char*gr,char*sth,char*pv);
#endif
