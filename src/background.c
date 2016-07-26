#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif

#include "background.h"
#include "sequence.h"

int getNucleotideFrequencyFromSequence(FILE *f, double *di, int order, int *nseq, int *lseq) {
  char *buffer=NULL;
  Sequence seq;
  int i=0, j;

  double ret=0;
 
 	allocSequence(&seq, *nseq, lseq);
 	getSequence(f, &seq);
  //memset(seq,0,(lmax+1)*sizeof(char));

  for (i=0; i<seq.nseq; i++) {
   	//if (seq.lseq[i]==0) break;
    for (j=0; j<seq.lseq[i]; j++) {
      if (j>=order) {
        di[getIndexFromAssignment(&seq.seq[i][j-order], order+1)]+=1.0;
        di[getIndexFromReverseAssignment(&seq.seq[i][j-order], order+1)]+=1.0;
        di[getIndexFromComplementaryAssignment(&seq.seq[i][j-order],order+1)]+=1.0;
        di[getIndexFromReverseComplementaryAssignment(&seq.seq[i][j-order],order+1)]+=1.0;
      }
    }
  }

  if (buffer) Free(buffer);
  destroySequence(&seq);

  return ret;
}

int getStationaryDistribution(double *trans, double *station, int order) {
  // use LAPACK to compute stationary distribution
  //
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
  	error("Memory-allocation in getStationaryDistribution failed");
	}
  tmpres=tmp1;
  tmpstart=tmp2;

  for (i=0; i<power(ALPHABETSIZE, order); i++) {
    tmpstart[i]=1.0/(double)power(ALPHABETSIZE, order);
  }

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

      //for (a=0; a<ALPHABETSIZE; a++) {
        tmpres[nextindex]+=tmpstart[previndex]*trans[j];
      //}
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

int getForwardTransition(double *di, double *trans, int order) {
  int i=0, j=0;
  double sum=0,ret=0;

  ret=sum;
  for (i=0; i<power(ALPHABETSIZE, order+1); i+=4) {
    sum=0.0;
    for (j=0; j< ALPHABETSIZE; j++) {
      sum+=di[i+j];
    }
    for (j=0; j< ALPHABETSIZE; j++) {
      trans[i+j]=di[i+j]/sum;
    }
  }
  return ret;
}

void writeBackground(FILE *f, double * station, double * trans, int order) {
  fwrite(&order, sizeof(int),1, f);
  if (order==0) {
    fwrite(station, sizeof(double), power(ALPHABETSIZE, order+1), f);
    fwrite(station, sizeof(double), power(ALPHABETSIZE, order+1), f);
  } else {
    fwrite(station, sizeof(double), power(ALPHABETSIZE, order), f);
    fwrite(trans, sizeof(double), power(ALPHABETSIZE, order+1), f);
  }
}

void deleteBackground(double * station, double *trans) {
  Free(station);
  if (trans!=NULL) {
    Free(trans);
  }
}

void readBackground (FILE *f, double **station, double **trans, int *order) {
  fread(order, sizeof(int),1, f);

  if (*order>0) {
    *station=Calloc(power(ALPHABETSIZE, *order), double);
    *trans=Calloc(power(ALPHABETSIZE, *order+1), double);
    if (*station==NULL||*trans==NULL) {
  	  error("Memory-allocation in readBackground failed");
	  }
    fread(*station, sizeof(double),power(ALPHABETSIZE, *order),f );
    fread(*trans, sizeof(double),power(ALPHABETSIZE, *order+1),f );
  } else {
    *station=Calloc(ALPHABETSIZE, double);
    *trans=Calloc(ALPHABETSIZE, double);
    if (*station==NULL||*trans==NULL) {
  	  error("Memory-allocation in readBackground failed");
	  }
    fread(*station, sizeof(double),ALPHABETSIZE,f );
    memcpy(*trans, *station, sizeof(double)*ALPHABETSIZE);
  }
}

void printBackground(double *stat, double *trans, int order) {
  int i;
	if (stat==NULL) return;

  if (order>0) {
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
      Rprintf("mu(i=%d)=%e\n", i, stat[i]);
    }
    for (i=0; i<power(ALPHABETSIZE, order+1); i++) {
      Rprintf("T(i=%d)=%e\n", i, trans[i]);
    }
	} else {

    for (i=0; i<power(ALPHABETSIZE, 1); i++) {
      Rprintf("mu(i=%d)=%e\n", i, stat[i]);
    }
    for (i=0; i<power(ALPHABETSIZE, 1); i++) {
      Rprintf("T(i=%d)=%e\n", i, trans[i]);
    }
	}
}

