#ifndef minmaxscore_h
#define minmaxscore_h
#include "matrix.h"

typedef struct {
    double dx; // is supposed to be the same as in struct ScoreMetaInfo
    int *maxforward;
    int *maxbackward;
    int *minforward;
    int *minbackward;
    int *intervalstart;
    int *intervalend;
    int Smax;
    int xmin;
    int order;
    int len;
} ExtremalScore;

int initExtremalScore(ExtremalScore *s, double, int length, int order);
int min(int a, int b);
int max(int a, int b);
int deleteExtremalScore(ExtremalScore *s);

void maxScoresPerPositionBack(DMatrix *theta, double *bg1, int *result, double *dx, int order);
void minScoresPerPositionBack(DMatrix *theta, double *bg1, int *result, double *dx, int order);
void maxScoresPerPositionForward(DMatrix *theta, double *bg1, double *trans, int *result, double *dx, int order);
void minScoresPerPositionForward(DMatrix *theta, double *bg1, double *trans, int *result, double *dx, int order);

void maxMotifScoreBack(DMatrix *,double *, int *, double *dx, int *ret, int order);
void minMotifScoreBack(DMatrix *,double *, int *, double *dx, int *ret, int order);
void loadMinMaxScores(DMatrix *theta, double *station, double *trans, ExtremalScore *loaded);
//void loadMinMaxScores(DMatrix *pwm, double *station, double *trans, ExtremalScore *e) {
void loadIntervalSize(ExtremalScore *escore, int *threshold);

//int maxScoreIntervalSize(ExtremalScore *escore);
int getScoreLowerBound(ExtremalScore *escore, int pos, int index);
int getScoreUpperBound(ExtremalScore *escore, int pos, int nucindex);
int getScoreLowerBoundUnconstrainted(ExtremalScore *escore, int pos, int nucindex);
int getScoreUpperBoundUnconstrainted(ExtremalScore *escore, int pos, int nucindex);

int maxScoreIntervalSize(ExtremalScore *e);
int getTotalScoreLowerBound(ExtremalScore *escore);
int getTotalScoreUpperBound(ExtremalScore *escore);
//int maxScorePerPositionWithThreshold(DMatrix *theta, double *bg1, ExtremalScore *result, int order, int threshold);
int * getScoreLowerBoundArray(ExtremalScore *escore, int pos);
int * getScoreUpperBoundArray(ExtremalScore *escore, int pos);
int * getLastScoreLowerBound(ExtremalScore *escore);
int getScoreLowerBoundPos(ExtremalScore *escore, int pos);
int getScoreUpperBoundPos(ExtremalScore *escore, int pos);
//int minScorePerPositionWithThreshold(DMatrix *theta, double *bg1, ExtremalScore *result, int order, int threshold);
int getMax(int *v, int N);
int getMin(int *v, int N);

void testmax(char *pwmfile, char *bgfile, char *output, char*);
void testmaxthreshold(char *pwmfile, char *bgfile, char *output, char *gran, char *threshold);

#endif
