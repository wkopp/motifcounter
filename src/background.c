#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif

#include "background.h"
#include "sequence.h"
#include "matrix.h"

void getNucleotideFrequencyFromSequence(char *seq, int slen,
        double *counts, int order) {
    int j;

    if (getSequenceLength(seq, slen)<0) {
        return;
    }

    for (j=order; j<slen; j++) {
        counts[getIndexFromAssignment(&seq[j-order], order+1)]+=1.0;
        counts[getIndexFromReverseAssignment(
                            &seq[j-order], order+1)]+=1.0;
        counts[getIndexFromComplementaryAssignment(
                            &seq[j-order], order+1)]+=1.0;
        counts[getIndexFromReverseComplementaryAssignment(
                            &seq[j-order],order+1)]+=1.0;
    }
}

int getStationaryDistribution(double *trans, double *station, int order) {
    int j, i, k;
    double *tmp1, *tmp2, *tmpres, *tmpstart;
    int ass [order+1];
    int nextindex, previndex;

    if (order <1) {
        error("no stationary distribution needed");
        return 1;
    }

    tmp1=Calloc(power(ALPHABETSIZE, order), double);
    tmp2=Calloc(power(ALPHABETSIZE, order), double);
    if (tmp1==NULL||tmp2==NULL) {
        error("Memory-allocation in getStationarydistribution failed");
    }
    tmpres=tmp1;
    tmpstart=tmp2;

    for (i=0; i<power(ALPHABETSIZE, order); i++) {
        tmpstart[i]=1.0/(double)power(ALPHABETSIZE, order);
    }

    // compute the stationary distribution using the power method
    //
    for (i=0; i<1000; i++) {
        for (j=0; j<power(ALPHABETSIZE, order+1); j++) {
            getAssignmentFromIndex(j, order+1, ass);
            nextindex=0;
            previndex=0;
            for (k=1; k<order+1; k++) {
                nextindex+=ass[k]*power(ALPHABETSIZE, order-k);
            }
            for (k=0; k<order; k++) {
                previndex+=ass[k]*power(ALPHABETSIZE, order-k-1);
            }

            tmpres[nextindex]+=tmpstart[previndex]*trans[j];
        }
        if (tmpres==tmp1) {
            tmpstart=tmpres;
            tmpres=tmp2;
        } else {
            tmpstart=tmpres;
            tmpres=tmp1;
        }
        memset(tmpres, 0, power(ALPHABETSIZE, order)* sizeof(double));
    }

    memmove(station,tmpstart, power(ALPHABETSIZE, order)* sizeof(double));
    Free(tmp1);
    Free(tmp2);
    return 0;
}

int getForwardTransition(double *counts, double *trans, int order) {
    int i=0, j=0;
    double sum=0,ret=0;

    ret=sum;
    for (i=0; i<power(ALPHABETSIZE, order+1); i+=4) {
        sum=0.0;
        for (j=0; j< ALPHABETSIZE; j++) {
            sum+=counts[i+j];
        }
        for (j=0; j< ALPHABETSIZE; j++) {
            trans[i+j]=counts[i+j]/sum;
        }
    }
    return ret;
}
