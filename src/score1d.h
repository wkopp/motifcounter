#ifndef score1d_h
#define score1d_h

#include "matrix.h"
#include "minmaxscore.h"
#include "scorefunctions.h"
#include "sequence.h"


typedef struct {
    double *y;
    double merged;
    int start, end;
} Score1d;

typedef struct {
    int mlen;
    ScoreMetaInfo meta;
    Score1d *ScoreBuffer1;
    Score1d *tmpScore;
    Score1d totalScore;
} MotifScore1d;

// initialize the score distribution structure
int initScoreDistribution1d (DMatrix *theta, double *bg1, MotifScore1d *result,
                             int order);

double getQuantileWithIndex1d(MotifScore1d *s, int qi);
int getQuantileIndex1d(Score1d *s, double pvalue);
double getProbWithIndex1d(MotifScore1d *s, int iquantile);
double getProb1d(MotifScore1d *s, double quantile);


int computeMarginalScoreDistribution1d(DMatrix *theta, double *bg1,
                                       double *bg0,
                                       MotifScore1d *marginalScore, int order);
int computeMarginalScoreDistribution1dBruteForce(DMatrix *pwm, double *trans,
        double *station, MotifScore1d *mscore, int xmin, int order);

int computeScoreDistribution1d(DMatrix *theta, double *bg1, double *bg0,
                               MotifScore1d *marginalScore,
                               ExtremalScore *,
                               int );

void scoredist(char *bgfile, char *pwmfile, char *oq, char *od, char *, char *,
               char *);
void scoredist_bruteforce(char *bgfile, char *pwmfile, char *output_quantile,
                          char *output_dist, char *pv, char *, char *gr);
void printscore(char *bgfile, char *pwmfile, char *gr);
void resetScore1d(Score1d *score, ScoreMetaInfo *meta);
double partionFunctionScore1d(Score1d *a, ScoreMetaInfo *meta);
double getProbability1d(Score1d *a, ScoreMetaInfo *meta);
void cutScoreRangeWithThreshold(MotifScore1d *mscore, ExtremalScore *tm,
                                int order);
#endif
