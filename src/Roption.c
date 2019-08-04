#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <string.h>

double Rsiglevel = 0.0, Rgran = 0.0;

int Rnucmap[26];
int Rrevcomp[4];

void Roption(double *siglevel, double *gran, int *ncores) {
    Rsiglevel = siglevel[0];
    Rgran = gran[0];

    memset(Rnucmap, -1, sizeof(Rnucmap));
    Rnucmap[0] = 0;
    Rnucmap['C'-'A'] = 1;
    Rnucmap['G'-'A'] = 2;
    Rnucmap['T'-'A'] = 3;

    Rrevcomp[0] = 3;
    Rrevcomp[1] = 2;
    Rrevcomp[2] = 1;
    Rrevcomp[3] = 0;

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
