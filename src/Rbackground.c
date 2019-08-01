#include <R.h>
#include <Rinternals.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"
#include "matrix.h"


SEXP Rcountfreq(SEXP rseqs, SEXP rorder) {
    int nseqs;
    const char *seq;
    int slen;
    int i;
    int *order = INTEGER(rorder);
    SEXP counts;
    double *xcounts;

    if (!isNewList(rseqs)) {
        Rprintf("Not a list. A list of sequence must be supplied.");
        return R_NilValue;
    }

    nseqs = length(rseqs);

    counts = PROTECT(allocVector(REALSXP, power(ALPHABETSIZE, order[0] + 1)));

    xcounts = REAL(counts);
    memset(xcounts, 0, (power(ALPHABETSIZE, order[0] + 1))*sizeof(double));

    for (i=0; i<nseqs; i++) {
      seq = CHAR(STRING_ELT(VECTOR_ELT(rseqs, i), 0));
      slen = strlen(seq);
      getNucleotideFrequencyFromSequence(seq, slen, xcounts, order[0]);
    }
    UNPROTECT(1);
    return counts;
}

void Rbgfromfreq(double *counts, double *station, double *trans, int *order) {
    int i;
    // check if all transition probabilities are greater than zero
    for (i = 0; i < power(ALPHABETSIZE, order[0] + 1); i++) {
        if (counts[i] <= 0) {
            error("All transition probabilities must be greater than zero:"
                  "Either reduce the order of the Markov model or use a DNA "
                  "sequence that is more heterogeneous");
        }
    }

    if (order[0] > 0) {
        getForwardTransition(counts, trans, order[0]);
        getStationaryDistribution(trans, station, order[0]);
    } else {
        getForwardTransition(counts, trans, order[0]);
        getForwardTransition(counts, station, order[0]);
    }
}
