#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sequence.h"
#include "matrix.h"
#include "scorefunctions.h"
#include "minmaxscore.h"

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
    IMatrix pwm;

    pfm.data = (double*)R_alloc((size_t)nrow[0] * ncol[0], sizeof(double));
    memset(pfm.data, 0, nrow[0] * ncol[0]*sizeof(double));

    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));


    pwm.nrow = pfm.nrow - order[0];
    pwm.ncol = power(ALPHABETSIZE, order[0] + 1);
    pwm.data = (int*)R_alloc((size_t) pwm.nrow * pwm.ncol, sizeof(int));
    memset(pwm.data, 0, pwm.nrow * pwm.ncol * sizeof(int));

    // from the PFM and background, produce a PWM
    getPositionWeights(station, trans, &pfm, &pwm, Rgran, order[0]);

    scores = PROTECT(allocVector(REALSXP, slen - pfm.nrow + 1));
    xscores = REAL(scores);

    scoreSequence(&pwm, seq, slen, xscores,
                  Rgran, order[0]);

    UNPROTECT(1);
    return scores;
}


SEXP Rhitsequence(SEXP rpfm_, SEXP rnrow, SEXP rncol, SEXP rseq,
                  SEXP rstation, SEXP rtrans, SEXP rorder, SEXP rthreshold) {
    double *pfm_ = REAL(rpfm_);
    double *station = REAL(rstation);
    double *trans = REAL(rtrans);
    int *nrow = INTEGER(rnrow);
    int *ncol = INTEGER(rncol);
    int *order = INTEGER(rorder);
    double *xhits;
    double *threshold = REAL(rthreshold);
    const char *seq;
    int slen;
    ExtremalScore escore;
  
    SEXP hits;
  
    seq = CHAR(STRING_ELT(rseq, 0));
    slen = strlen(seq);
  
    if (getSequenceLength(seq, slen) < ncol[0]) {
       return R_NilValue;
    }
  
    DMatrix pfm;
    IMatrix pwm;
  
    pfm.data = (double*)R_alloc((size_t)nrow[0] * ncol[0], sizeof(double));
    memset(pfm.data, 0, nrow[0] * ncol[0]*sizeof(double));
  
    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));
  
    pwm.nrow = pfm.nrow - order[0];
    pwm.ncol = power(ALPHABETSIZE, order[0] + 1);
    pwm.data = (int*)R_alloc((size_t) pwm.nrow * pwm.ncol, sizeof(int));
    memset(pwm.data, 0, pwm.nrow * pwm.ncol * sizeof(int));
  
    // from the PFM and background, produce a PWM
    getPositionWeights(station, trans, &pfm, &pwm, Rgran, order[0]);
  
    hits = PROTECT(allocVector(REALSXP, slen - pfm.nrow + 1));
    xhits = REAL(hits);
    memset(xhits, 0, (slen - pfm.nrow +1)*sizeof(double));
    
    // fill up max min scores
    // allocate memory
    initExtremalScore(&escore, Rgran, pfm.nrow, order[0]);
    
    // load min max values per position
    loadMinMaxScores(&pfm, station, trans, &escore);
    
    hitSequence(&pwm, seq, slen, xhits,
                Rgran, order[0], threshold[0], &escore);
  
    UNPROTECT(1);
    return hits;
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

    // otherwise, compute the histogram
    scoreHistogram(station, trans,
                   &pfm, seq, slen, xdist, Rgran, mins, order[0]);


    UNPROTECT(2);
    return dist;
}
