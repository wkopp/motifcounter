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

extern int Rorder;
extern int RorderForSampling;
extern double *Rstation, *Rtrans;
extern double *RstationForSampling, *RtransForSampling;
extern double Rsiglevel, Rgran;

void RgenRndSeq(char **seq, int *len) {
    if (!RstationForSampling) {
        error("load background properly");
    }
    generateRandomSequence(RstationForSampling, RtransForSampling, seq[0],
        len[0], RorderForSampling);
}

