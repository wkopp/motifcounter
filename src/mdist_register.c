#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>


void Rcountfreq(char **seq, int *slen, double *counts, int *order);
void Rbgfromfreq(double *counts, double *station, double *trans, int *order);
static R_NativePrimitiveArgType countfreq_t[] = {
   STRSXP, INTSXP,REALSXP, INTSXP
};
static R_NativePrimitiveArgType bgfromfreq_t[] = {
   REALSXP, REALSXP,REALSXP, INTSXP
};
void RPosteriorProbability(double *alpha, double *beta,
          double *beta3p, double *beta5p,
            double *hitdistribution, int *sseqlen,
              int *smaxhits, int *snos,int *motiflen,
              int *singlestranded);

static R_NativePrimitiveArgType combinatorial_t[] = {
   REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
   INTSXP,INTSXP,INTSXP,INTSXP,INTSXP
};

void Roverlap(double *data,int*nrow,int*ncol,double *alpha, double *beta,
        double *beta3p, double *beta5p, double *gamma,
    double *,double *,int*);
void RoverlapSingleStranded(double *data,int*nrow,int*ncol,
        double *alpha, double *beta,
        double *beta3p, double *beta5p, double *gamma,
    double *,double *,int*);
static R_NativePrimitiveArgType overlap_t[] = {
        REALSXP,INTSXP,INTSXP,REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
        REALSXP,REALSXP,INTSXP
};

void RcompoundpoissonPape_useGamma(double *gamma,
          double *hitdistribution, int *nseq,
          int *lseq, int * mhit, int *mclump, int *motiflen);

static R_NativePrimitiveArgType cp_gamma_t[] = {
        REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,INTSXP,INTSXP
};
void Rcompoundpoisson_useBeta(double *alpha, double *beta,
          double *beta3p, double *beta5p,
            double *hitdistribution, int *nseq, int *lseq, int * mhit,
            int *mclump, int *motiflen, int *sstrand);
static R_NativePrimitiveArgType cp_beta_t[] = {
        REALSXP, REALSXP, REALSXP,REALSXP,
        REALSXP,INTSXP, INTSXP, INTSXP,INTSXP,INTSXP, INTSXP
};
void RgenRndSeq(char **num_seqs, int *len,
double *,double *,int*);
static R_NativePrimitiveArgType num_seqs_t[] = {
        STRSXP,INTSXP,REALSXP,REALSXP,INTSXP
};

void Roption(double *siglevel, double *gran, int *ncores);
static R_NativePrimitiveArgType option_t[] = {
    REALSXP,REALSXP,INTSXP
};
void Rfsiglevel(double *siglevel);
static R_NativePrimitiveArgType siglevel_t[] = {
    REALSXP
};
void Rscorerange(double *,int*,int*,int *scorerange,
double *,double *,int*);
static R_NativePrimitiveArgType scorerange_t[] = {
    REALSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,INTSXP
};
void Rscoredist(double *,int*,int*, double *score, double *prob,
double *,double *,int*);
void Rscoredist_bf(double *,int*,int*, double *score, double *prob,
double *,double *,int*);
static R_NativePrimitiveArgType scoredist_t[] = {
    REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP
};
void Rscoresequence(double *,int*,int*, char **,
    double *fscore,double *rscore, int *len,
double *,double *,int*);
static R_NativePrimitiveArgType scoreseq_t[] = {
    REALSXP,INTSXP,INTSXP,STRSXP,REALSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP
};
void RscoreHistogram(double *,int*,int*,char **, int *,
            double *scores, double *freq,
        double *,double *,int*);
static R_NativePrimitiveArgType scorehist_t[] = {
    REALSXP,INTSXP,INTSXP,STRSXP,INTSXP,REALSXP,REALSXP
    ,REALSXP,REALSXP,INTSXP
};

static R_CMethodDef cMethods[] = {
    {"motifcounter_countfreq", (DL_FUNC) &Rcountfreq, 4, countfreq_t},
    {"motifcounter_bgfromfreq", (DL_FUNC) &Rbgfromfreq, 4, bgfromfreq_t},
    {"motifcounter_combinatorialDist",
            (DL_FUNC) &RPosteriorProbability, 10, combinatorial_t},
    {"motifcounter_overlapSingleStranded",
            (DL_FUNC) &RoverlapSingleStranded, 11, overlap_t},
    {"motifcounter_overlap", (DL_FUNC) &Roverlap, 11, overlap_t},
    {"motifcounter_compoundPoisson_useBeta",
            (DL_FUNC) &Rcompoundpoisson_useBeta, 11, cp_beta_t},
    {"motifcounter_compoundPoissonPape_useGamma",
            (DL_FUNC) &RcompoundpoissonPape_useGamma, 7, cp_gamma_t},
    {"motifcounter_generateRndSeq", (DL_FUNC) &RgenRndSeq, 5, num_seqs_t},
    {"motifcounter_option", (DL_FUNC) &Roption, 3, option_t},
    {"motifcounter_siglevel", (DL_FUNC) &Rfsiglevel, 1, siglevel_t},
    {"motifcounter_scorerange", (DL_FUNC) &Rscorerange, 7, scorerange_t},
    {"motifcounter_scoredist", (DL_FUNC) &Rscoredist, 8, scoredist_t},
    {"motifcounter_scoresequence", (DL_FUNC) &Rscoresequence, 10, scoreseq_t},
    {"motifcounter_scoredist_bf", (DL_FUNC) &Rscoredist_bf, 8, scoredist_t},
    {"motifcounter_scorehistogram",
            (DL_FUNC) &RscoreHistogram, 10, scorehist_t},
    {NULL, NULL, 0}
};



void R_init_motifcounter(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

void R_unload_motifcounter(DllInfo *info) {}
