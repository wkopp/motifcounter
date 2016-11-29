#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "background.h"
#include "sequence.h"
#include "scorefunctions.h"
#include "minmaxscore.h"
#include "score1d.h"

char sampleNucleotide(double *prob) {
    double ip=unif_rand();
    double p;
    double cumprob=0.0;
    char n;
    int i;

    p=(double)(ip);
    for (i=0; i<ALPHABETSIZE; i++) {
        cumprob+=prob[i];
        if(p<=cumprob) {
            break;
        }
    }
    if (i==4)
        n= getNuc(i-1);
    else
        n= getNuc(i);
    return n;
}

void sampleInitialNucleotide(double *prob, char *seq, int order) {
    double ip=unif_rand();
    double p;
    double cumprob=0.0;
    int i, ass[order];

    p=(double)(ip);
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
        cumprob+=prob[i];
        if(p<=cumprob) {
            break;
        }
    }
    getAssignmentFromIndex(i, order, ass);
    for(i=0;i<order; i++) {
        seq[i]=getNuc(ass[i]);
    }
}

void generateRandomSequence(double *station, double *trans, char *seq,
        int seqlen, int order) {
    int i, ind;

    sampleInitialNucleotide(station, seq, order);
    for (i=order; i<seqlen; i++) {
        ind=getIndexFromAssignment(&seq[i-order], order);
        seq[i]=sampleNucleotide(&trans[ind*ALPHABETSIZE]);
    }
}

