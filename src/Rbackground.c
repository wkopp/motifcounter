#include <R.h>
#include <Rinternals.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"

double *Rstation=NULL, *Rtrans=NULL;
int Rorder;
void RdestroyBackground();

void Rmakebg(char **infasta, int *order, int *nseq, int *lseq) {
    FILE *f;
    double *count;
    int i;

    RdestroyBackground();
    f =fopen(infasta[0],"r");
    if (f==NULL) {
        error("%s not found!",infasta[0]);
        return;
    }

    if (order[0]>0) {
        Rstation=Calloc(power(ALPHABETSIZE,order[0]),double);
        count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        Rtrans=Calloc(power(ALPHABETSIZE, order[0]+1),double);
        if (Rstation==NULL || count==NULL || Rtrans==NULL) {
            error("Memory allocation in Rmakebg failed");
        }

        getNucleotideFrequencyFromSequence(f,count, order[0], nseq, lseq);

        getForwardTransition(count, Rtrans, order[0]);
        getStationaryDistribution(Rtrans, Rstation, order[0]);

    } else {
        Rstation=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        Rtrans=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        if (Rstation==NULL || count==NULL || Rtrans==NULL) {
            error("Memory allocation in Rmakebg failed");
        }
        getNucleotideFrequencyFromSequence(f,count, order[0], nseq,lseq);
        getForwardTransition(count, Rstation, order[0]);
        getForwardTransition(count, Rtrans, order[0]);

    }
    // check if all transition probabilities are greater than zero
    for (i=0; i<power(ALPHABETSIZE,order[0]+1);i++) {
        if (count[i]<=0) {
            RdestroyBackground();
            fclose(f);
            Free(count);
            error("All transition probabilities must be greater than zero:"
                    "Either reduce the order of the Markov model or use a DNA "
                    "sequence that is more heterogeneous");
            return;
        }
    }


    fclose(f);
    Rorder=order[0];

    Free(count);
}


void RgetBackground(double *station, double *trans) {
    int i;
    if (Rstation==NULL) return;
    for (i=0; i<4; i++) {
        Rprintf("%e\n", Rstation[i]);
    }
    for (i=0; i<4; i++) {
        Rprintf("%e\n", Rtrans[i]);
    }
}

void RprintBackground() {
    printBackground(Rstation,Rtrans,Rorder);
}

void RdestroyBackground() {
    if(Rstation) Free(Rstation);
    if(Rtrans) Free(Rtrans);
    Rstation=NULL;
    Rtrans=NULL;
}

SEXP fetchStationBackground() {
    int i,or_;
    SEXP station;
    if (!Rstation) return R_NilValue;

    or_=(Rorder==0) ? 1 : Rorder;

    PROTECT(station=allocVector(REALSXP,power(ALPHABETSIZE,or_)));
    for (i=0; i<power(ALPHABETSIZE,or_); i++) {
        REAL(station)[i]=Rstation[i];
    }
    UNPROTECT(1);
    return station;
}

SEXP fetchTransBackground() {
    int i,or_;
    SEXP trans;
    if (!Rtrans) return R_NilValue;

    or_=Rorder+1;

    PROTECT(trans=allocVector(REALSXP,power(ALPHABETSIZE,or_)));
    for (i=0; i<power(ALPHABETSIZE,or_); i++) {
        REAL(trans)[i]=Rtrans[i];
    }
    UNPROTECT(1);
    return trans;
}
