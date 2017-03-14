#ifndef overlap_h
#define overlap_h
#include "score1d.h"

void computeDeltas(double *delta, double *deltap,
                   double *beta, double *beta3p, double *beta5p, int mlen);
void computeDeltasSingleStranded(double *delta,
                                 double *beta,  int mlen);
void computeBetas(double *beta, double *beta3p, double *beta5p,
                  double *gamma, int mlen, double eps);
void computeBetasSingleStranded(double *beta,
                                double *gamma, int mlen, double eps);
void overlap(char *bgfile, char *pwmfile, char *output_quantile,
             char *output_dist, char *pv, char *gr, char *method);
#endif
