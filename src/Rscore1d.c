#include <R.h>
#include <Rinternals.h>
#include "score1d.h"


extern double Rsiglevel, Rgran;

SEXP Rscoredist(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                SEXP rstation, SEXP rtrans, SEXP rorder) {
    MotifScore1d null;
    ExtremalScore fescore;
    double dx;
    double *xdist;
    double *pfm_ = REAL(rpfm_);
    double *station = REAL(rstation);
    double *trans = REAL(rtrans);
    int *nrow = INTEGER(rnrow);
    int *ncol = INTEGER(rncol);
    int *order = INTEGER(rorder);
    SEXP dist;
    int i, intervalsize, maxs, mins;
    DMatrix pfm;

    pfm.data = Calloc(nrow[0] * ncol[0], double);

    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));

    dx = Rgran;

    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins = getTotalScoreLowerBound(&fescore);
    maxs = getTotalScoreUpperBound(&fescore);

    intervalsize = maxs - mins;

    initScoreMetaInfo(mins,
                      maxs, intervalsize, dx, &null.meta);
    null.meta.prob = &ProbBg;
    null.meta.probinit = &ProbinitBg;
    initScoreDistribution1d(&pfm, trans, &null, order[0]);

    computeScoreDistribution1d(&pfm, trans,  station,
                               &null, &fescore, order[0]);

    dist = PROTECT(allocVector(REALSXP, null.meta.xmax - null.meta.xmin + 1));
    xdist = REAL(dist);
    for (i = 0; i < null.meta.xmax - null.meta.xmin + 1; i++) {
        xdist[i] = null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
    UNPROTECT(1);
    return dist;
}

SEXP Rscoredist_bf(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                   SEXP rstation, SEXP rtrans, SEXP rorder) {
    int i;
    MotifScore1d null;
    ExtremalScore fescore;
    double dx;
    double *xdist;
    double *pfm_ = REAL(rpfm_);
    double *station = REAL(rstation);
    double *trans = REAL(rtrans);
    int *nrow = INTEGER(rnrow);
    int *ncol = INTEGER(rncol);
    int *order = INTEGER(rorder);
    int mins, maxs;
    int intervalsize;
    SEXP dist;
    DMatrix pfm;

    pfm.data = Calloc(nrow[0] * ncol[0], double);

    // Rcol and c-col are swapped
    pfm.ncol = nrow[0];
    pfm.nrow = ncol[0];
    memcpy(pfm.data, pfm_, nrow[0]*ncol[0]*sizeof(double));

    dx = Rgran;

    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins = getTotalScoreLowerBound(&fescore);
    maxs = getTotalScoreUpperBound(&fescore);

    intervalsize = maxs - mins;

    initScoreMetaInfo(mins,
                      maxs, intervalsize, dx, &null.meta);

    null.meta.prob = &ProbBg;
    null.meta.probinit = &ProbinitBg;
    initScoreDistribution1d(&pfm, trans, &null, order[0]);

    computeMarginalScoreDistribution1dBruteForce(&pfm, trans,
            station, &null, null.meta.xmin, order[0]);

    dist = PROTECT(allocVector(REALSXP, null.meta.xmax - null.meta.xmin + 1));
    xdist = REAL(dist);
    for (i = 0; i < null.meta.xmax - null.meta.xmin + 1; i++) {
        xdist[i] = null.totalScore.y[i];
    }
    deleteExtremalScore(&fescore);
    deleteScoreDistribution1d(&null, order[0]);

    Free(pfm.data);
    UNPROTECT(1);
    return dist;
}
