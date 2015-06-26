#ifndef markovchain_h
#define markovchain_h
void markovchain(double *dist, 
  double *alpha, double *beta, double *beta3p, double *beta5p, int slen);
double minmc(int n, double *par, void *extra);
void dmc(int n, double *par, double *gradient, void *ex);
void removeDist();
#endif
