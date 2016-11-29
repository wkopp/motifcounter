#include <stdio.h>
#include <R.h>
#include "sequence.h"
#include "matrix.h"

int isNucleotide(char c) {
    int r=0;
    switch (c) {
        case 'a':
        case 'A':
        case 'c':
        case 'C':
        case 'g':
        case 'G':
        case 't':
        case 'T':
            r=1;
            break;
        case 'N':
        case 'n':
            r=-1;
            break;
    }
    return r;
}

int getNucIndex(char c) {
    int r=-1;
    switch (c) {
        case 'a':
        case 'A':
            r=0;
            break;
        case 'c':
        case 'C':
            r=1;
            break;
        case 'g':
        case 'G':
            r=2;
            break;
        case 't':
        case 'T':
            r=3;
            break;
    }
    return r;
}

int getComplementNucIndex(char c) {
    int r=-1;
    switch (c) {
        case 'a':
        case 'A':
            r=3;
            break;
        case 'c':
        case 'C':
            r=2;
            break;
        case 'g':
        case 'G':
            r=1;
            break;
        case 't':
        case 'T':
            r=0;
            break;
        case 'n':
        case 'N':
            r=-1;
            break;
    }
    return r;
}

char getNuc(int i) {
    char r='n';
    switch (i) {
        case 0:
            r='a';
            break;
        case 1:
            r='c';
            break;
        case 2:
            r='g';
            break;
        case 3:
            r='t';
            break;
    }
    return r;
}


int getSequence(FILE *f, Sequence *s) {
    int writeseq=0, writeheader=0, i=0;
    int iseq=-1, ipos=0;
    int lmax=0;
    char *buffer;
    double ret=0;

    for (i=0; i< s->nseq; i++) {
        if (lmax<s->lseq[i]) {
            lmax=s->lseq[i];
        }
    }
    buffer=Calloc(lmax+1, char);
    if(buffer==NULL) {
        error("Memory-allocation in getSequence failed");
    }

    while(fgets(buffer, sizeof(buffer), f)!=NULL) {
        for (i=0; i<strlen(buffer); i++) {
            if (buffer[i]=='>') {
                iseq++;
                if (s->lseq[iseq]>0) {
                    writeheader=1;
                }  else {
                    writeheader=0;
                }
                writeseq=0;
            }
            if (writeheader==1 && buffer[i]=='\n') {
                writeheader=0; writeseq=1; ipos=0; break;
            }

            if (writeseq==1 && buffer[i]=='\n') break;
            if (writeseq==1 && isNucleotide(buffer[i])==1) {
                s->seq[iseq][ipos++]=buffer[i];
            }
        }
    }
    Free(buffer);

    return ret;
}

void allocSequence(Sequence *seq, int nseq, int *lseq) {
    int i;
    seq->seq=Calloc(nseq, char*);
    seq->lseq=Calloc(nseq, int);
    if (seq->seq==NULL||seq->lseq==NULL) {
        error("Memory-allocation in allocSequence failed");
    }
    for (i=0; i<nseq; i++) {
  	seq->seq[i]=Calloc(lseq[i]+1,char);
        if (seq->seq[i]==NULL) {
            error("Memory-allocation in allocSequence failed");
        }
        seq->lseq[i]=lseq[i];
    }
    seq->nseq=nseq;
}

void destroySequence(Sequence *seq) {
    int i;
    for (i=0; i<seq->nseq; i++) {
  	Free(seq->seq[i]);
    }
    Free(seq->lseq);
    Free(seq->seq);
}


void getAssignmentFromIndex(int index, int length, int *ret) {
    int rest=index;
    int i;

    for (i=0; i<length; i++) {
        ret[i]=rest/power(ALPHABETSIZE,length-1-i);
        rest-=ret[i]*power(ALPHABETSIZE,length-1-i);
    }
    return;
}

int getIndexFromAssignment(char *ass, int length) {
    int index=0;
    int i;

    for (i=0; i<length; i++) {
        index+=getNucIndex(ass[i])*power(ALPHABETSIZE,length-1-i);
    }
    return index;
}

int getIndexFromReverseAssignment(char *ass, int length) {
    int index=0;
    int i;

    for (i=0; i<length; i++) {
        index+=getNucIndex(ass[length-1-i])*power(ALPHABETSIZE,length-1-i);
    }
    return index;
}

int getIndexFromComplementaryAssignment(char *ass, int length) {
    int index=0;
    int i;

    for (i=0; i<length; i++) {
        index+=getComplementNucIndex(ass[i]) *power(ALPHABETSIZE,length-1-i);
    }
    return index;
}

int getIndexFromReverseComplementaryAssignment(char *ass, int length) {
    int index=0;
    int i;

    for (i=0; i<length; i++) {
        index+=getComplementNucIndex(ass[length-1-i]) *
            power(ALPHABETSIZE,length-i-1);
    }
    return index;
}
