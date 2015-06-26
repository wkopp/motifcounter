#ifndef countdist_h
#define countdist_h

void deleteCountDistribution(double *cp);
double* initCountDistribution(int maxhit);
void writeCountDistribution(FILE *f, double *cp, int maxhit);
#endif
