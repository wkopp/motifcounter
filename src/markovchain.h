#ifndef markovchain_h
#define markovchain_h
void markovchain(double *dist,
                 double *alpha, double *beta, double *beta3p,
                 double *beta5p, int slen, int motiflen);
double minmc(int n, double *par, void *extra);
void dmc(int n, double *par, double *gradient, void *ex);
void removeDist();

typedef struct {
    double alpha;
    double *beta;
    double *beta3p;
    double *beta5p;
    int len;
    int motiflen;
} CGParams;
#endif
