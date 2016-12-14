#include <R.h>
#include <Rinternals.h>
#include "minmaxscore.h"


extern double Rgran;

SEXP Rscorerange(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                    SEXP rstation, SEXP rtrans, SEXP rorder) {
    ExtremalScore fescore;
    int mins, maxs;
    double dx;
    double *xscores;
    double *pfm_=REAL(rpfm_);
    double *station=REAL(rstation);
    double *trans=REAL(rtrans);
    int *nrow=INTEGER(rnrow);
    int *ncol=INTEGER(rncol);
    int *order=INTEGER(rorder);
    SEXP scores;
    int i;
    DMatrix pfm;

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    if (pfm.data==NULL) {
        error("Rscorerange: memory allocation failed");
    }
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));

    dx=Rgran;
    initExtremalScore(&fescore, dx, pfm.nrow, order[0]);

    loadMinMaxScores(&pfm, station, trans, &fescore);
    loadIntervalSize(&fescore, NULL);

    mins=getTotalScoreLowerBound(&fescore);
    maxs=getTotalScoreUpperBound(&fescore);

    deleteExtremalScore(&fescore);

    Free(pfm.data);
    
    scores=PROTECT(allocVector(REALSXP, maxs-mins + 1));
    xscores=REAL(scores);
    for (i=0; i<maxs-mins + 1; i++) {
      xscores[i]=(double)(mins+i)*dx;
    }
    
    UNPROTECT(1);
    return scores;
}
