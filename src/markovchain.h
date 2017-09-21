#ifndef markovchain_h
#define markovchain_h


double getOptimalTauMCDS(double *alpha, double *beta, double *beta3p,
    double *beta5p, int *motiflen);
double getOptimalTauMCSS(double *alpha, double *beta, int *motiflen);

typedef struct {
    double alpha;
    double *beta;
    double *beta3p;
    double *beta5p;
    int len;
    int motiflen;
    double * dist;
} CGParams;
#endif
