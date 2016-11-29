#include <R.h>
#include "sequence.h"
#include "matrix.h"
#include "scorefunctions.h"

extern double Rgran;
extern double *Rstation, *Rtrans;
extern int Rorder;

void RnumSeqs(char ** fastafile, int *numofseqs) {
    FILE *f;
    char buffer[1024*16];
    numofseqs[0]=0;
    f =fopen(fastafile[0],"r");
    if (f==NULL) {
        error("IO-Error in RnumSeqs during opening of fasta file");
    }

    while(fgets(buffer, sizeof(buffer), f)!=NULL) {
        if (buffer[0]=='>') {
            numofseqs[0]++;
        }
    }
    if (ferror(f)) {
        error("IO-Error in RnumSeqs");
    }
    fclose(f);
}

void RlenSeqs(char ** fastafile, int *numofseqs, int * lseq) {
    FILE *f;
    char buffer[1024*16];
    int i, iseq=0;
    int writeheader=0, writeseq=0;
    f =fopen(fastafile[0],"r");
    if (f==NULL) {
        error("Error in RnumSeqs opening fasta file");
    }

    while(fgets(buffer, sizeof(buffer), f)!=NULL) {
        for (i=0; i<strlen(buffer); i++) {
            if (buffer[i]=='>') {
                lseq[iseq]=0;
                iseq++;
                writeheader=1;
                writeseq=0;
            }
            if (writeseq==1 && buffer[i]=='\n') break;
            if (writeseq==1 && isNucleotide(buffer[i])==1) lseq[iseq-1]++;
            if (writeseq==1 && isNucleotide(buffer[i])<0) {
                lseq[iseq-1]=0;
                //warning("Sequence number %d contains 'n' or"
                //"  'N' and is discarded.",iseq);
                writeseq=0;
                break;
            }
            if (writeheader==1 && buffer[i]=='\n') {
                writeheader=0;
                writeseq=1;
                break;
            }
        }
    }
    if (iseq!= numofseqs[0]) {
        error("RlenSeqs: Number of sequences does not match!");
    }
    if (ferror(f)) {
        error("IO-Error in RnumSeqs");
    }
    fclose(f);
}


void Rscoresequence(double *pfm_, int *nrow, int *ncol, char **seq,
    double *fscores, double *rscores, int *slen) {
    int i, n;

    double dx;
    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }
    if (Rstation==NULL || Rtrans==NULL) {
        error("Background model uninitialized! "
                "Use readBackground()");
        return;
    }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    cpfm.data=Calloc(nrow[0]*ncol[0],double);
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    cpfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    cpfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));
    for (i=1; i<=nrow[0]*ncol[0];i++) {
        cpfm.data[i-1]=pfm.data[nrow[0]*ncol[0]-i];
    }

    dx=Rgran;
    scoreSequence(Rstation, Rtrans,
        &pfm, seq[0], slen[0], fscores,
        Rgran, Rorder);
    scoreSequence(Rstation, Rtrans,
        &cpfm, seq[0], slen[0], rscores,
        Rgran, Rorder);

    Free(pfm.data);
    Free(cpfm.data);
}


void RscoreHistogram(double *pfm_, int *nrow, int *ncol, char **seq, int *slen,
    double *scorebins,  double *frequency) {
    int i, n;
    ExtremalScore fescore, rescore;
    int fmins,fmaxs, rmins,rmaxs;
    int mins, maxs, noscores;
    DMatrix pfm, cpfm;

    if (Rgran==0.0) {
        error("call mdistOption  first");
        return;
    }
    if (Rstation==NULL || Rtrans==NULL) {
        error("Background model uninitialized! "
                "Use readBackgroundForSampling()");
        return;
    }

    pfm.data=Calloc(nrow[0]*ncol[0],double);
    cpfm.data=Calloc(nrow[0]*ncol[0],double);
    // Rcol and c-col are swapped
    pfm.ncol=nrow[0];
    cpfm.ncol=nrow[0];
    pfm.nrow=ncol[0];
    cpfm.nrow=ncol[0];
    memcpy(pfm.data,pfm_,nrow[0]*ncol[0]*sizeof(double));
    for (i=1; i<=nrow[0]*ncol[0];i++) {
        cpfm.data[i-1]=pfm.data[nrow[0]*ncol[0]-i];
    }

    initExtremalScore(&fescore, Rgran, pfm.nrow, Rorder);
    initExtremalScore(&rescore, Rgran, cpfm.nrow, Rorder);

    loadMinMaxScores(&pfm, Rstation, Rtrans, &fescore);
    loadMinMaxScores(&cpfm, Rstation, Rtrans, &rescore);
    loadIntervalSize(&fescore, NULL);
    loadIntervalSize(&rescore, NULL);

    fmins=getTotalScoreLowerBound(&fescore);
    rmins=getTotalScoreLowerBound(&rescore);
    fmaxs=getTotalScoreUpperBound(&fescore);
    rmaxs=getTotalScoreUpperBound(&rescore);
    maxs=(fmaxs>rmaxs) ? fmaxs : rmaxs;
    mins=(fmins>rmins) ? fmins : rmins;

    // if the sequence contains any N's, do not process the scores
    noscores=0;
    for (i=0; i<slen[0];i++) {
        if (getNucIndex(seq[0][i])<0) {
            noscores=1;
            break;
        }
    }
    if (noscores==0) {
        scoreHistogram(Rstation, Rtrans,
            &pfm, seq[0], slen[0], frequency, Rgran, mins,Rorder);
    }
    for (i=0; i<maxs-mins+1; i++) {
        scorebins[i]= (double)(mins+i)*Rgran;
    }
    deleteExtremalScore(&fescore);
    deleteExtremalScore(&rescore);

    Free(pfm.data);
    Free(cpfm.data);
}
