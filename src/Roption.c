#ifdef _OPENMP
#include <omp.h>
#endif

double Rsiglevel=0.0, Rgran=0.0;

void Roption(double *siglevel, double *gran, int *ncores) {
    Rsiglevel=siglevel[0];
    Rgran=gran[0];
#ifdef _OPENMP
    omp_set_num_threads(*ncores);
#endif
}
