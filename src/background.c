#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif

#include "background.h"
#include "sequence.h"

#ifdef WK
int getSeqen(FILE *f) {
  char dig;
  int writeseq=0, writeheader=1, i=0;

  while((dig=fgetc(f))!=EOF) {
//Rprintf("%c",dig);
    //if (dig=='>') {
    //  writeheader=1;
    //}
    if (writeseq==1 && dig=='\n') continue;
    if (writeheader==1 && dig=='\n') {
      writeheader=0;
      writeseq=1;
      continue;
    }
    if (writeseq==1) {
      i++;
    }
  }
  //Rprintf("len=%d\n", i);
  return i;
}
#endif

int getNucleotideFrequencyFromSequence(FILE *f, double *di, int order, int *nseq, int *lseq) {
  char *buffer=NULL, *seq=NULL;
  int writeseq=0, writeheader=0, i=0, j;
  int iseq;
  int lmax=-1;

  double ret=0;
 
  for (i=0; i< *nseq; i++) {
    if (lmax<lseq[i]) {
      lmax=lseq[i];
    }
  }
  buffer=Calloc(lmax+1, char);
  seq=Calloc(lmax+1, char);
  //memset(seq,0,(lmax+1)*sizeof(char));

  while(fgets(buffer, sizeof(buffer), f)!=NULL) {
    for (i=0; i<strlen(buffer); i++) {
      if (buffer[i]=='>') {
        memset(seq,0,(lmax+1)*sizeof(char));
        //lseq[iseq]=0;
        iseq=0;
        //iseq++;
        writeheader=1;
        writeseq=0;
      }
      if (writeheader==1 && buffer[i]=='\n') { writeheader=0; writeseq=1; iseq=0; break; }
      if (writeseq==1 && buffer[i]=='\n') break;
      if (writeseq==1 && isNucleotide(buffer[i])==1) seq[iseq++]=buffer[i];
      if (writeseq==1 && isNucleotide(buffer[i])<0) {
        //lseq[iseq-1]=0;
        warning("Sequence number %d contains 'n' or 'N' and is discarded.",iseq);
        memset(seq,0,(lmax+1)*sizeof(char));
        writeseq=0;
        break;
      }
      if (buffer[i]=='>' && writeseq==1) {
      
        for (j=0; j<strlen(seq); j++) {
          if (j>=order && isNucleotide(seq[j])==1) {
            di[getIndexFromAssignment(&seq[j-order], order+1)]+=1.0;
            di[getIndexFromReverseAssignment(&seq[j-order], order+1)]+=1.0;
            di[getIndexFromComplementaryAssignment(&seq[j-order],order+1)]+=1.0;
            di[getIndexFromReverseComplementaryAssignment(&seq[j-order],order+1)]+=1.0;
          }
        }
      }
    }
  }
  if (feof(f) && writeseq==1) {
    for (j=0; j<strlen(seq); j++) {
      if (j>=order && isNucleotide(seq[j])==1) {
        di[getIndexFromAssignment(&seq[j-order], order+1)]+=1.0;
        di[getIndexFromReverseAssignment(&seq[j-order], order+1)]+=1.0;
        di[getIndexFromComplementaryAssignment(&seq[j-order],order+1)]+=1.0;
        di[getIndexFromReverseComplementaryAssignment(&seq[j-order],order+1)]+=1.0;
      }
    }
	}

  if (buffer) Free(buffer);
  if (seq) Free(seq);

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
  #ifdef IN_R
    error("no stationary distribution needed");
  #else
    fprintf(stderr, "no stationary distribution needed");
  #endif
    return 1;
  }

#ifdef IN_R
  tmp1=Calloc(power(ALPHABETSIZE, order), double);
  tmp2=Calloc(power(ALPHABETSIZE, order), double);
  #else
  tmp1=calloc(power(ALPHABETSIZE, order), sizeof(double));
  tmp2=calloc(power(ALPHABETSIZE, order), sizeof(double));
  #endif
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
  #ifdef IN_R
  Free(tmp1);
  Free(tmp2);
  #else
  free(tmp1);
  free(tmp2);
  #endif
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
#ifdef IN_R
  Free(station);
  if (trans!=NULL) {
    Free(trans);
  }
  #else
  free(station);
  if (trans!=NULL) {
    free(trans);
  }
  #endif
}

void readBackground (FILE *f, double **station, double **trans, int *order) {
  fread(order, sizeof(int),1, f);

  if (*order>0) {
  #ifdef IN_R
    *station=Calloc(power(ALPHABETSIZE, *order), double);
    *trans=Calloc(power(ALPHABETSIZE, *order+1), double);
    #else
    *station=calloc(power(ALPHABETSIZE, *order), sizeof(double));
    *trans=calloc(power(ALPHABETSIZE, *order+1), sizeof(double));
    #endif
    fread(*station, sizeof(double),power(ALPHABETSIZE, *order),f );
    fread(*trans, sizeof(double),power(ALPHABETSIZE, *order+1),f );
  } else {
  #ifdef IN_R
    *station=Calloc(ALPHABETSIZE, double);
    *trans=Calloc(ALPHABETSIZE, double);
    #else
    *station=calloc(ALPHABETSIZE, sizeof(double));
    *trans=calloc(ALPHABETSIZE, sizeof(double));
    #endif
    fread(*station, sizeof(double),ALPHABETSIZE,f );
    memcpy(*trans, *station, sizeof(double)*ALPHABETSIZE);
  }
}

#ifndef IN_R
int makebg(char *infasta, char *outmodel, char *outseq, char *order) {
  FILE *f;
  double *stationary,*count, *trans=NULL;
  int iorder;
  #ifdef DEBUG
  int i,j;
  #endif
  //int seqlen;
  //Sequence seq;

    if (order==NULL) {
      iorder=1;
    } else {
      iorder=atol(order);
    }
    f =fopen(infasta,"r");
    if (f==NULL) {
      fprintf(stderr, "error during file opening");
      return 1;
    }


    if (iorder>0) {
    #ifdef IN_R
      stationary=Calloc(power(ALPHABETSIZE,iorder),double);
      count=Calloc(power(ALPHABETSIZE,iorder+1),double);
      trans=Calloc(power(ALPHABETSIZE, iorder+1),double);
      #else
      stationary=calloc(power(ALPHABETSIZE,iorder),sizeof(double));
      count=calloc(power(ALPHABETSIZE,iorder+1),sizeof(double));
      trans=calloc(power(ALPHABETSIZE, iorder+1),sizeof(double));
      #endif

      getNucleotideFrequencyFromSequence(f,count, iorder);

      getForwardTransition(count, trans, iorder);
      getStationaryDistribution(trans, stationary, iorder);

    } else {

    #ifdef IN_R
      stationary=Calloc(power(ALPHABETSIZE,iorder+1),double);
      count=Calloc(power(ALPHABETSIZE,iorder+1),double);
      #else
      stationary=calloc(power(ALPHABETSIZE,iorder+1),sizeof(double));
      count=calloc(power(ALPHABETSIZE,iorder+1),sizeof(double));
      #endif
      getNucleotideFrequencyFromSequence(f,count, iorder);
      getForwardTransition(count, stationary, iorder);
    }
    fclose(f);

#ifdef DEBUG
    printf(stdout,"Counts:\n");
    for (j=0; j < power(ALPHABETSIZE, iorder+1); j++) {
      printf("%1.5e ", count[j]);
      printf("\n");
    }
    if (iorder>0) {
      printf("Stationary distribution\n");
      for (i=0; i < power(ALPHABETSIZE, iorder); i++) {
        printf( "%1.7e ", stationary[i]);
        printf("\n");
      }
      printf("Transition-Probabilities:\n");
      for (j=0; j < power(ALPHABETSIZE, iorder+1); j++) {
        printf( "%1.7e ", trans[j]);
        printf("\n");
      }
    } else {
      printf("Stationary distribution\n");
      for (i=0; i < power(ALPHABETSIZE, iorder+1); i++) {
        printf( "%1.7e ", stationary[i]);
        printf("\n");
      }
    }
    #endif
    #ifdef IN_R
    Free(count);
    #else
    free(count);
    #endif
    f =fopen(outmodel,"ab");
    if (f==NULL) {
      deleteBackground(stationary, trans);
      fprintf(stderr, "error during file opening");
      return 1;
    }
    //Write the length of the DHS into the file before
    // the background parameters
    //fwrite(&seqlen, sizeof(int),1,f);
    writeBackground(f, stationary, trans, iorder);
    deleteBackground(stationary, trans);
    fclose(f);

    return 0;
}
#endif
