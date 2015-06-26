#include <stdio.h>
#include <stdlib.h>
#ifdef IN_R
#include <R.h>
#endif

#ifndef IN_R
double *initCountDistribution(int maxhit) {
#ifdef IN_R
  return Calloc(maxhit+1, double);
  #else
  return calloc(maxhit+1, sizeof(double));
  #endif
}

void deleteCountDistribution(double *cp) {
#ifdef IN_R
  Free(cp);
  #else
  free(cp);
  #endif
}

void writeCountDistribution(FILE *f, double *cp, int maxhit) {
  int i;
  for (i=0;i<maxhit+1;i++) {
    fprintf(f, "%d\t%1.4e\n",i,cp[i]);
  }
}
#endif
