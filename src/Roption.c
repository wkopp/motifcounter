#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>

double Rsiglevel = 0.0, Rgran = 0.0;

void Roption(double *siglevel, double *gran, int *ncores) {
    Rsiglevel = siglevel[0];
    Rgran = gran[0];
#ifdef _OPENMP
    omp_set_num_threads(*ncores);
#endif
}
void Rfsiglevel(double *siglevel) {
    if (Rsiglevel <= 0.0) {
        error("Uninitialized alpha. Call motifcounterOption first");
    }
    siglevel[0] = Rsiglevel;
}
