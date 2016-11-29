#include <R.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"

double *RstationForSampling=NULL, *RtransForSampling=NULL;
int RorderForSampling;
void RdestroyBackgroundForSampling();

void RmakebgForSampling(char **infasta, int *order, int *nseq, int *lseq) {
    FILE *f;
    double *count;
    int i;

    RdestroyBackgroundForSampling();
    f =fopen(infasta[0],"r");
    if (f==NULL) {
        error("%s not found!",infasta[0]);
        return;
    }


    if (order[0]>0) {
        RstationForSampling=Calloc(power(ALPHABETSIZE,order[0]),double);
        count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        RtransForSampling=Calloc(power(ALPHABETSIZE, order[0]+1),double);
        if (RstationForSampling==NULL || count==NULL ||
                    RtransForSampling==NULL) {
            error("Memory allocation in RmakebgForSampling failed");
        }

        getNucleotideFrequencyFromSequence(f,count, order[0], nseq, lseq);
        getForwardTransition(count, RtransForSampling, order[0]);
        getStationaryDistribution(RtransForSampling,
                    RstationForSampling, order[0]);

    } else {
        RstationForSampling=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        RtransForSampling=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
        if (RstationForSampling==NULL || count==NULL ||
                    RtransForSampling==NULL) {
            error("Memory allocation in RmakebgForSampling failed");
        }
        getNucleotideFrequencyFromSequence(f,count, order[0], nseq, lseq);
        getForwardTransition(count, RstationForSampling, order[0]);
        getForwardTransition(count, RtransForSampling, order[0]);
    }
    // check if all transition probabilities are greater than zero
    for (i=0; i<power(ALPHABETSIZE,order[0]+1);i++) {
        if (count[i]<=0) {
            RdestroyBackgroundForSampling();
            fclose(f);
            Free(count);
            error("All transition probabilities must be greater than zero:"
                    "Either reduce the order of the Markov model or use a DNA "
                    "sequence that is more heterogeneous");
        }
    }

    fclose(f);
    RorderForSampling=order[0];

    Free(count);
}

void RgetOrderForSampling(int *o) {
    if (!RstationForSampling) {
        error("Load background for sampling first!");
    }
    o[0]=RorderForSampling;
}
void RgetBackgroundForSampling(double *station, double *trans) {
    int i;
    for (i=0; i<4; i++) {
        Rprintf("%e\n", RstationForSampling[i]);
    }
    for (i=0; i<4; i++) {
        Rprintf("%e\n", RtransForSampling[i]);
    }
}

void RprintBackgroundForSampling() {
    printBackground(RstationForSampling,RtransForSampling,RorderForSampling);
}

void RdestroyBackgroundForSampling() {
    if(RstationForSampling) Free(RstationForSampling);
    if(RtransForSampling) Free(RtransForSampling);
    RstationForSampling=NULL;
    RtransForSampling=NULL;
}
