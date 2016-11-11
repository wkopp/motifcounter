#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fetchStationBackground();
SEXP fetchTransBackground();
SEXP fetchMotif();

static R_CallMethodDef callMethods[]  = {
       {"mdist_fetchStationBackground", (DL_FUNC) &fetchStationBackground, 0},
       {"mdist_fetchTransBackground", (DL_FUNC) &fetchTransBackground, 0},
       {"mdist_fetchMotif", (DL_FUNC) &fetchMotif, 0},
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
              int *smaxhits, int *snos);

static R_NativePrimitiveArgType combinatorial_t[] = {
   REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,INTSXP,INTSXP,INTSXP,INTSXP
};

void Roverlap(double *alpha, double *beta, double *beta3p, double *beta5p, double *gamma);
void RoverlapSingleStranded(double *alpha, double *beta, double *beta3p, double *beta5p,
                        double *gamma);
static R_NativePrimitiveArgType overlap_t[] = {
        REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

void RcompoundpoissonPape_useGamma(double *gamma,
          double *hitdistribution, int *nseq, int *lseq, int * mhit, int *mclump);

static R_NativePrimitiveArgType cp_gamma_t[] = {
        REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,INTSXP
};
void Rcompoundpoisson_useBeta(double *alpha, double *beta,
          double *beta3p, double *beta5p,
            double *hitdistribution, int *nseq, int *lseq, int * mhit, 
            int *mclump, int *sstrand);
static R_NativePrimitiveArgType cp_beta_t[] = {
        REALSXP, REALSXP, REALSXP,REALSXP,REALSXP,INTSXP, INTSXP, INTSXP,INTSXP,INTSXP
};
void Rloadmotif(double *data, int *nrow, int *ncol);
static R_NativePrimitiveArgType loadmotif_t[] = {
        REALSXP, INTSXP, INTSXP
};
void Rmotiffromfile(char **fmotif, double *pseudocount);
static R_NativePrimitiveArgType motiffromfile_t[] = {
        STRSXP, REALSXP
};
void Rdestroymotif();
void Rmotiflength(int *mlen);
static R_NativePrimitiveArgType motif_len_t[] = {
        INTSXP
};
void RnumSeqs(char ** fastafile, int *numofseqs);
static R_NativePrimitiveArgType num_seqs_t[] = {
        STRSXP,INTSXP
};
void RlenSeqs(char ** fastafile, int *numofseqs, int * lseq);
static R_NativePrimitiveArgType seq_len_t[] = {
        STRSXP,INTSXP,INTSXP
};
void RnumberOfHits(char **inputfile, 
        int *numofhits, int *nseq, int *lseq, int *singlestranded);
static R_NativePrimitiveArgType num_hits_t[] = {
    STRSXP,INTSXP,INTSXP,INTSXP,INTSXP
};
void Roption(double *siglevel, double *gran, int *ncores);
static R_NativePrimitiveArgType option_t[] = {
    REALSXP,REALSXP,INTSXP
};
void Rscorerange(int *scorerange);
static R_NativePrimitiveArgType scorerange_t[] = {
    INTSXP
};
void Rscoredist( double *score, double *prob);
void Rscoredist_bf( double *score, double *prob);
static R_NativePrimitiveArgType scoredist_t[] = {
    REALSXP,REALSXP
};
void RsimulateScores(double *scores, double *distribution, int *slen,
          int *perm);
static R_NativePrimitiveArgType simscoredist_t[] = {
    REALSXP,REALSXP,INTSXP,INTSXP
};

void RsimulateCountDistribution( double *distribution, int* perm,
           int *nseq, int *lseq, int *mxhit, int *singlestranded);
static R_NativePrimitiveArgType simcountdist_t[] = {
    REALSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP
};

static R_CMethodDef cMethods[] = {
    {"mdist_makebg", (DL_FUNC) &Rmakebg, 4, makebg_t},
    {"mdist_printBackground", (DL_FUNC) &RprintBackground, 0, NULL},
    {"mdist_deleteBackground", (DL_FUNC) &RdestroyBackground, 0, NULL},
    {"mdist_makebgForSampling", (DL_FUNC) &RmakebgForSampling, 4, makebg_t},
    {"mdist_printBackgroundForSampling", (DL_FUNC) &RprintBackgroundForSampling, 0, NULL},
    {"mdist_deleteBackgroundForSampling", (DL_FUNC) &RdestroyBackgroundForSampling, 0, NULL},
    {"mdist_combinatorialDist", (DL_FUNC) &RPosteriorProbability, 9, combinatorial_t},
    {"mdist_overlapSingleStranded", (DL_FUNC) &RoverlapSingleStranded, 5, overlap_t},
    {"mdist_overlap", (DL_FUNC) &Roverlap, 5, overlap_t},
    {"mdist_compoundPoisson_useBeta", (DL_FUNC) &Rcompoundpoisson_useBeta, 10, cp_beta_t},
    {"mdist_compoundPoissonPape_useGamma", (DL_FUNC) &RcompoundpoissonPape_useGamma, 6, cp_gamma_t},
    {"mdist_loadmotif", (DL_FUNC) &Rloadmotif, 3, loadmotif_t},
    {"mdist_motiffromfile", (DL_FUNC) &Rmotiffromfile, 2, motiffromfile_t},
    {"mdist_deleteMotif", (DL_FUNC) &Rdestroymotif, 0, NULL},
    {"mdist_motiflength", (DL_FUNC) &Rmotiflength, 1, motif_len_t},
    {"mdist_numSeqs", (DL_FUNC) &RnumSeqs, 2, num_seqs_t},
    {"mdist_lenSeqs", (DL_FUNC) &RlenSeqs, 3, seq_len_t},
    {"mdist_numberOfHits", (DL_FUNC) &RnumberOfHits, 5, num_hits_t},
    {"mdist_option", (DL_FUNC) &Roption, 3, option_t},
    {"mdist_scorerange", (DL_FUNC) &Rscorerange, 1, scorerange_t},
    {"mdist_scoredist", (DL_FUNC) &Rscoredist, 2, scoredist_t},
    {"mdist_scoredist_bf", (DL_FUNC) &Rscoredist_bf, 2, scoredist_t},
    {"mdist_simulateScores", (DL_FUNC) &RsimulateScores, 4, simscoredist_t},
    {"mdist_simulateCountDistribution", (DL_FUNC) &RsimulateCountDistribution, 6, simcountdist_t},
    {NULL, NULL, 0}
};



void R_init_mdist(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}


