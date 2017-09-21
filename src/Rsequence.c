#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sequence.h"
#include "matrix.h"
#include "scorefunctions.h"

extern double Rgran;

SEXP Rslen(SEXP rseq) {
    const char *seq;
    int slen;

    seq = CHAR(STRING_ELT(rseq, 0));
    slen = strlen(seq);

    slen = getSequenceLength(seq, slen);

    if (slen < 0) {
        return ScalarInteger(0);
    } else {
        return ScalarInteger(slen);
    }
}

SEXP Rscoresequence(SEXP rpfm_, SEXP rnrow, SEXP rncol, SEXP rseq,
                    SEXP rstation, SEXP rtrans, SEXP rorder) {
    double *pfm_ = REAL(rpfm_);
    double *station = REAL(rstation);
    double *trans = REAL(rtrans);
    int *nrow = INTEGER(rnrow);
    int *ncol = INTEGER(rncol);
    int *order = INTEGER(rorder);
    double *xscores;
    const char *seq;
    int slen;

    SEXP scores;

    seq = CHAR(STRING_ELT(rseq, 0));
    slen = strlen(seq);

    if (getSequenceLength(seq, slen) < ncol[0]) {
        return R_NilValue;
    }

    DMatrix pfm;

    pfm.data = (double*)R_alloc((size_t)nrow[0] * ncol[0], sizeof(double));
    memset(pfm.data, 0, nrow[0] * ncol[0]*sizeof(double));

    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));

    scores = PROTECT(allocVector(REALSXP, slen - pfm.nrow + 1));
    xscores = REAL(scores);

    scoreSequence(station, trans,
                  &pfm, seq, slen, xscores,
                  Rgran, order[0]);

    UNPROTECT(1);
    return scores;
}


SEXP RscoreHistogram(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                     SEXP rseq, SEXP rstation, SEXP rtrans, SEXP rorder) {
    int i;
    ExtremalScore fescore;
    int mins, maxs, noscores;
    DMatrix pfm;
    double *pfm_ = REAL(rpfm_);
    double *station = REAL(rstation);
    double *trans = REAL(rtrans);
    const char *seq;
    int slen;
    int *nrow = INTEGER(rnrow);
    int *ncol = INTEGER(rncol);
    int *order = INTEGER(rorder);
    SEXP dist;
    double *xdist;

    PROTECT(rseq = AS_CHARACTER(rseq));
    seq = CHAR(STRING_ELT(rseq, 0));
    slen = strlen(seq);

    pfm.data = (double*)R_alloc((size_t)nrow[0] * ncol[0], sizeof(double));
    memset(pfm.data, 0, nrow[0] * ncol[0]*sizeof(double));

    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));

    initExtremalScore(&fescore, Rgran, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins = getTotalScoreLowerBound(&fescore);
    maxs = getTotalScoreUpperBound(&fescore);

    dist = PROTECT(allocVector(REALSXP, maxs - mins + 1));
    xdist = REAL(dist);
    for (i = 0; i < maxs - mins + 1; i++) xdist[i] = 0.0;


    noscores = 0;
    if (getSequenceLength(seq, slen) < 0) {
        noscores = 1;
    }

    // if there is a non-nucleotide, just return the empty histogram
    if (noscores == 0) {
        // otherwise, compute the histogram
        scoreHistogram(station, trans,
                       &pfm, seq, slen, xdist, Rgran, mins, order[0]);
    }

    UNPROTECT(2);
    return dist;
}
