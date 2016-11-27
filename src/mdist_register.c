#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fetchStationBackground();
SEXP fetchTransBackground();
//SEXP fetchMotif();

static R_CallMethodDef callMethods[]  = {
    {"motifcounter_fetchStationBackground", 
        (DL_FUNC) &fetchStationBackground, 0},
    {"motifcounter_fetchTransBackground", (DL_FUNC) &fetchTransBackground, 0},
    //{"motifcounter_fetchMotif", (DL_FUNC) &fetchMotif, 0},
    {NULL, NULL, 0}
};



void Rmakebg(char **infasta, int *order, int *nseq, int *lseq);
static R_NativePrimitiveArgType makebg_t[] = {
        STRSXP, INTSXP, INTSXP, INTSXP
};
void RprintBackground();
void RdestroyBackground();
void RmakebgForSampling(char **infasta, int *order, int *nseq, int *lseq);
void RprintBackgroundForSampling();
void RdestroyBackgroundForSampling();

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
        double *beta3p, double *beta5p, double *gamma);
void RoverlapSingleStranded(double *data,int*nrow,int*ncol,
        double *alpha, double *beta, 
        double *beta3p, double *beta5p, double *gamma);
static R_NativePrimitiveArgType overlap_t[] = {
        REALSXP,INTSXP,INTSXP,REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
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
void RnumSeqs(char ** fastafile, int *numofseqs);
static R_NativePrimitiveArgType num_seqs_t[] = {
        STRSXP,INTSXP
};
void RlenSeqs(char ** fastafile, int *numofseqs, int * lseq);
static R_NativePrimitiveArgType seq_len_t[] = {
        STRSXP,INTSXP,INTSXP
};
void RnumberOfHits(double *pfm, int*nrow,int*ncol,char **inputfile, 
        int *numofhits, int *nseq, int *lseq, int *singlestranded);
static R_NativePrimitiveArgType num_hits_t[] = {
    REALSXP,INTSXP,INTSXP,STRSXP,INTSXP,INTSXP,INTSXP,INTSXP
};
void Roption(double *siglevel, double *gran, int *ncores);
static R_NativePrimitiveArgType option_t[] = {
    REALSXP,REALSXP,INTSXP
};
void Rscorerange(double *,int*,int*,int *scorerange);
static R_NativePrimitiveArgType scorerange_t[] = {
    REALSXP,INTSXP,INTSXP,INTSXP
};
void Rscoredist(double *,int*,int*, double *score, double *prob);
void Rscoredist_bf(double *,int*,int*, double *score, double *prob);
static R_NativePrimitiveArgType scoredist_t[] = {
    REALSXP,INTSXP,INTSXP,REALSXP,REALSXP
};
void RsimulateScores(double *,int*,int*,double *scores, 
            double *distribution, int *slen,
            int *perm);
static R_NativePrimitiveArgType simscoredist_t[] = {
    REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,INTSXP,INTSXP
};

void RsimulateCountDistribution(double*,int*,int*, 
        double *distribution, int* perm,
        int *nseq, int *lseq, int *mxhit, int *singlestranded);
static R_NativePrimitiveArgType simcountdist_t[] = {
    REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP
};

static R_CMethodDef cMethods[] = {
    {"motifcounter_makebg", (DL_FUNC) &Rmakebg, 4, makebg_t},
    {"motifcounter_printBackground", (DL_FUNC) &RprintBackground, 0, NULL},
    {"motifcounter_deleteBackground", (DL_FUNC) &RdestroyBackground, 0, NULL},
    {"motifcounter_makebgForSampling", 
        (DL_FUNC) &RmakebgForSampling, 4, makebg_t},
    {"motifcounter_printBackgroundForSampling", 
            (DL_FUNC) &RprintBackgroundForSampling, 0, NULL},
    {"motifcounter_deleteBackgroundForSampling", 
            (DL_FUNC) &RdestroyBackgroundForSampling, 0, NULL},
    {"motifcounter_combinatorialDist", 
            (DL_FUNC) &RPosteriorProbability, 10, combinatorial_t},
    {"motifcounter_overlapSingleStranded", 
            (DL_FUNC) &RoverlapSingleStranded, 8, overlap_t},
    {"motifcounter_overlap", (DL_FUNC) &Roverlap, 8, overlap_t},
    {"motifcounter_compoundPoisson_useBeta", 
            (DL_FUNC) &Rcompoundpoisson_useBeta, 11, cp_beta_t},
    {"motifcounter_compoundPoissonPape_useGamma", 
            (DL_FUNC) &RcompoundpoissonPape_useGamma, 7, cp_gamma_t},
    {"motifcounter_numSeqs", (DL_FUNC) &RnumSeqs, 2, num_seqs_t},
    {"motifcounter_lenSeqs", (DL_FUNC) &RlenSeqs, 3, seq_len_t},
    {"motifcounter_numberOfHits", (DL_FUNC) &RnumberOfHits, 8, num_hits_t},
    {"motifcounter_option", (DL_FUNC) &Roption, 3, option_t},
    {"motifcounter_scorerange", (DL_FUNC) &Rscorerange, 4, scorerange_t},
    {"motifcounter_scoredist", (DL_FUNC) &Rscoredist, 5, scoredist_t},
    {"motifcounter_scoredist_bf", (DL_FUNC) &Rscoredist_bf, 5, scoredist_t},
    {"motifcounter_simulateScores", 
            (DL_FUNC) &RsimulateScores, 7, simscoredist_t},
    {"motifcounter_simulateCountDistribution", 
            (DL_FUNC) &RsimulateCountDistribution, 9, simcountdist_t},
    {NULL, NULL, 0}
};



void R_init_motifcounter(DllInfo *info) {
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

void R_unload_motifcounter(DllInfo *info) {
    RdestroyBackground();
    RdestroyBackgroundForSampling();
}
