#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>



void Rbgfromfreq(double *counts, double *station, double *trans, int *order);

static R_NativePrimitiveArgType bgfromfreq_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP
};
void RPosteriorProbability(double *alpha, double *beta,
                           double *beta3p, double *beta5p,
                           double *hitdistribution, int *sseqlen,
                           int *smaxhits, int *snos, int *motiflen,
                           int *singlestranded);

static R_NativePrimitiveArgType combinatorial_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

void Roverlap(double *data, int *nrow, int *ncol, double *alpha, double *beta,
              double *beta3p, double *beta5p, double *gamma,
              double *, double *, int *);
void RoverlapSingleStranded(double *data, int *nrow, int *ncol,
                            double *alpha, double *beta,
                            double *beta3p, double *beta5p, double *gamma,
                            double *, double *, int *);
static R_NativePrimitiveArgType overlap_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, INTSXP
};

void RcompoundpoissonPape_useGamma(double *gamma,
                                   double *hitdistribution, int *nseq,
                                   int *lseq, int *mhit, int *mclump, int *motiflen);

static R_NativePrimitiveArgType cp_gamma_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
void Rcompoundpoisson_useBeta(double *alpha, double *beta,
                              double *beta3p, double *beta5p,
                              double *hitdistribution, int *nseq, int *lseq, int *mhit,
                              int *mclump, int *motiflen, int *sstrand);
static R_NativePrimitiveArgType cp_beta_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
void RgenRndSeq(char **num_seqs, int *len,
                double *, double *, int *);
static R_NativePrimitiveArgType num_seqs_t[] = {
    STRSXP, INTSXP, REALSXP, REALSXP, INTSXP
};

void Roption(double *siglevel, double *gran, int *ncores);
static R_NativePrimitiveArgType option_t[] = {
    REALSXP, REALSXP, INTSXP
};
void Rfsiglevel(double *siglevel);
static R_NativePrimitiveArgType siglevel_t[] = {
    REALSXP
};

void RclumpsizeBeta(double *beta, double *beta3p, double *beta5p,
                        double *dist, int *maxclump, int *motiflen);

static R_NativePrimitiveArgType clumpsize_beta_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

void RclumpsizeGamma(double *gamma, double *dist, int *maxclump, int *motiflen);
static R_NativePrimitiveArgType clumpsize_gamma_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP
};

void markovchain_ss(double *dist, double *alpha, double *beta,
        int *steps, int *motiflen);
static R_NativePrimitiveArgType mcss_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

void markovchain(double *dist, double *alpha, double *beta,
        double *beta3p, double *beta5p, int *steps, int *motiflen);
static R_NativePrimitiveArgType mcds_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};


static R_CMethodDef cMethods[] = {

    {"motifcounter_bgfromfreq", (DL_FUNC) &Rbgfromfreq, 4, bgfromfreq_t},
    {"motifcounter_markovmodel_ss", (DL_FUNC) &markovchain_ss, 5, mcss_t},
    {"motifcounter_markovmodel_ds", (DL_FUNC) &markovchain, 7, mcds_t},

    {
        "motifcounter_combinatorialDist",
        (DL_FUNC) &RPosteriorProbability, 10, combinatorial_t
    },
    {
        "motifcounter_overlapSingleStranded",
        (DL_FUNC) &RoverlapSingleStranded, 11, overlap_t
    },
    {"motifcounter_overlap", (DL_FUNC) &Roverlap, 11, overlap_t},
    {
        "motifcounter_compoundPoisson_useBeta",
        (DL_FUNC) &Rcompoundpoisson_useBeta, 11, cp_beta_t
    },
    {
        "motifcounter_compoundPoissonPape_useGamma",
        (DL_FUNC) &RcompoundpoissonPape_useGamma, 7, cp_gamma_t
    },
    {"motifcounter_generateRndSeq", (DL_FUNC) &RgenRndSeq, 5, num_seqs_t},
    {"motifcounter_option", (DL_FUNC) &Roption, 3, option_t},
    {"motifcounter_siglevel", (DL_FUNC) &Rfsiglevel, 1, siglevel_t},
    {"motifcounter_clumpsize_kopp", (DL_FUNC) &RclumpsizeBeta, 6, clumpsize_beta_t},
    {"motifcounter_clumpsize_pape", (DL_FUNC) &RclumpsizeGamma, 4, clumpsize_gamma_t},
    {NULL, NULL, 0}
};

SEXP Rcountfreq(SEXP seqs, SEXP order);
SEXP RgetPositionWeights(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                    SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP Rscoresequence(SEXP rpfm_, SEXP rnrow, SEXP rncol, SEXP rseq,
                    SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP Rhitsequence(SEXP rpfm_, SEXP rnrow, SEXP rncol, SEXP rseq,
                  SEXP rstation, SEXP rtrans, SEXP rorder, SEXP rthreshold);
SEXP Rscorerange(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                 SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP Rscoredist(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP Rscoredist_bf(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                   SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP RscoreHistogram(SEXP rpfm_, SEXP rnrow, SEXP rncol,
                     SEXP rseq, SEXP rstation, SEXP rtrans, SEXP rorder);
SEXP Rslen(SEXP rseq);

SEXP mcss_check_optimal(SEXP alpha_, SEXP beta_, SEXP motiflen_);
SEXP mcds_check_optimal(SEXP alpha_, SEXP beta_, SEXP beta3p_,
    SEXP beta5p_, SEXP motiflen_);

SEXP Rmatchcount(SEXP rpfm_, SEXP rnrow, SEXP rncol, SEXP seqlist,
                 SEXP rstation, SEXP rtrans, SEXP rorder, SEXP rthreshold,
                 SEXP rignore_ns);

static R_CallMethodDef callMethods[]  = {
    {"motifcounter_countfreq", (DL_FUNC) &Rcountfreq, 2},
    {"motifcounter_slen", (DL_FUNC) &Rslen, 1},
    {"motifcounter_getpositionweights", (DL_FUNC) &RgetPositionWeights, 6},
    {"motifcounter_scoresequence", (DL_FUNC) &Rscoresequence, 7},
    {"motifcounter_hitsequence", (DL_FUNC) &Rhitsequence, 8},
    {"motifcounter_matchcount", (DL_FUNC) &Rmatchcount, 9},
    {"motifcounter_scorerange", (DL_FUNC) &Rscorerange, 6},
    {"motifcounter_scoredist", (DL_FUNC) &Rscoredist, 6},
    {"motifcounter_scoredist_bf", (DL_FUNC) &Rscoredist_bf, 6},
    {"motifcounter_scorehistogram", (DL_FUNC) &RscoreHistogram, 7},
    {"motifcounter_mcss_check_optimal", (DL_FUNC) &mcss_check_optimal, 3},
    {"motifcounter_mcds_check_optimal", (DL_FUNC) &mcds_check_optimal, 5},
    {NULL, NULL, 0}
};


void R_init_motifcounter(DllInfo *info) {
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

void R_unload_motifcounter(DllInfo *info) {}
