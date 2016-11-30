#include <time.h>
#include <R.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h"
#include "simulate.h"
#include "minmaxscore.h"
#include "background.h"
#include "score1d.h"


void RgenRndSeq(char **seq, int *len, double *station,
    double *trans, int *order) {

    generateRandomSequence(station, trans, seq[0],
        len[0], order[0]);
}
