#include <omp.h>

double Rsiglevel=0.0, Rgran=0.0;

void Roption(double *siglevel, double *gran, int *ncores) {
  Rsiglevel=siglevel[0];
  Rgran=gran[0];
  omp_set_num_threads(*ncores);
}
