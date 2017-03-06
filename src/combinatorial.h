#ifndef posteriorcount_h
#define posteriorcount_h

#define EPSILON 1e-30

typedef struct {
    int seqlen;
    int mlen;
    int maxhits;
    double ** *value;
    double *beta;
    double *beta3p;
    double *beta5p;
    double *delta;
    double *deltap;
    double alpha;
    double omega;
    double probzerohits;
    double **extentionMatrix;
    double *K;
    double *B;
} PosteriorCount;

double dpower(double v, int exp);
void testPosteriorProbability(char *bg, char *fpwm, char *output,
                              char *sseqlen,
                              char *smaxhits,  char *spv, char *sth, char *snos, char *gran);
int allocPosteriorProbability(PosteriorCount *p, int seqlen, int mlen,
                              int maxhits);
void deletePosteriorProbability(PosteriorCount *p);
void initPosteriorProbability(PosteriorCount *p, double alpha, double **beta,
                              double **beta3p, double **beta5p, double **delta, double **deltap);
void computePosteriorProbability(PosteriorCount *prob);
void finishPosteriorProbability(PosteriorCount *prob, double *final,
                                int nhits);
void multipleShortSequenceProbability(double *simple, double *aggregated,
                                      int maxsimplehits, int maxagghits, int numofseqs);
#endif
