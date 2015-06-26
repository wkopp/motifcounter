#ifndef compoundpoisson_h
#define compoundpoisson_h
#include <stdio.h>

double computePoissonParameter(int seqlen, int mlen, int maxclump, double alpha,double *theta);
double *initTheta(int maxclump);
void deleteTheta(double *theta);
void computeTheta(int maxclump, double *theta, double *delta, int mlen);
void computePapeTheta(int maxclump, double *theta, double *delta, double *gamma, int mlen);
void computeExtentionFactorsPape(double *xi, double *gamma, int mlen);
void computeExtentionFactorsKopp(double *xi,
  double *delta, double *deltap, double *beta, 
  double *beta3p, double *beta5p, int mlen);
void computeInitialClump(double *theta, double *gamma, int mlen);
void computeInitialClumpKopp(double *theta, double *gamma, 
  double *delta, double *deltap, int mlen);
double *initCompoundPoissonDistribution(int maxhit);
void computeCompoundPoissonDistribution(double lambda, 
  int maxhit, int maxclump, double *theta, double *cp);
void compoundpoisson_correctedpape(char *bgfile, char *pwmfile,
   char *output_dist, char *slen, char * mhit, char *mclump, char *pv, char *sth, char *gr);
void compoundpoisson_kopp(char *bgfile, char *pwmfile,
   char *output_dist, char *slen, char * mhit, char *mclump, char *pv, char *sth, char *gr);

#endif
