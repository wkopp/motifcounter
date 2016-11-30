#include <R.h>
#include <Rinternals.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"
#include "matrix.h"


void Rcountfreq(char **seq, int *slen, double *counts, int *order) {
    getNucleotideFrequencyFromSequence(seq[0],slen[0],counts,order[0]);
}

void Rbgfromfreq(double *counts, double *station, double *trans, int *order) {
    int i;
    // check if all transition probabilities are greater than zero
    for (i=0; i<power(ALPHABETSIZE,order[0]+1);i++) {
        if (counts[i]<=0) {
            error("All transition probabilities must be greater than zero:"
                    "Either reduce the order of the Markov model or use a DNA "
                    "sequence that is more heterogeneous");
        }
    }

    if (order[0]>0) {
        getForwardTransition(counts, trans, order[0]);
        getStationaryDistribution(trans, station, order[0]);
    } else {
        getForwardTransition(counts, trans, order[0]);
        getForwardTransition(counts, station, order[0]);
    }
}
