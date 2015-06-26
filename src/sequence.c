#include <stdio.h>
#include "sequence.h"
#include "matrix.h"

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
    default:
      r='n';
  }
  return r;
}


int getIndex(int i) {
  return i;
}

int getComplementFromIndex(int i) {
  int r=0;
  switch (i) {
    case 0:
      r=3;
      break;
    case 1:
      r=2;
      break;
    case 2:
      r=1;
      break;
    case 3:
      r=0;
      break;
    default:
      r=-1;
  }
  return r;
}

void writeSequenceBinary(FILE *f, Sequence *s, int dhs) {
  int i=0;
  //int dhs=DHSSIZE;
  fwrite(&s->numseq,sizeof(int),1,f);
  fwrite(&dhs,sizeof(int),1,f);
  for (i=0; i<s->numseq;i++) {
    fwrite(&s->seq[i][0], sizeof(char),dhs, f);
  }
}

int getSequenceFromBinary(FILE *f,Sequence *s) {
  int ret=0,dhs, i;
  ret=fread(&s->numseq,sizeof(int),1,f);
  if (ret<1) return -1;
  ret=fread(&dhs,sizeof(int),1,f);
  if (ret<1) return -1;
  //m->data=malloc(sizeof(double)*m->nrow*m->ncol);
  for (i=0; i<s->numseq; i++) {
    ret=fread(&s->seq[i][0],sizeof(char),dhs,f);
    if (ret<dhs) return -1;
  }
  return 0;
}

int getSequence(FILE *f, Sequence *s) {
  char dig;
  int writeseq=0, writeheader=0, i=0;
  double ret=0;
  while((dig=fgetc(f))!=EOF) {
    if (dig=='>') {
      s->numseq++;
      writeheader=1;
      i=0;
    }
    if (writeseq==1 && dig=='\n') writeseq=0;
    if (writeheader==1 && dig=='\n') {
      writeheader=0;
      writeseq=1;
      i=0;
      continue;
    }
    if (writeseq==1) {
      s->seq[s->numseq-1][i++]=dig;
    }
  }

  return ret;
}

int skipAssignment(char *ass, int len) {
  int i;
  for (i=0; i<len; i++) {
    if (ass[i]=='n' || ass[i]=='N') return 1;
  }
  return 0;
}

void getAssignmentFromIndex(int index, int length, int *ret) {
  int rest=index;
  int i;

  for (i=0; i<length; i++) {
   // if (rest <powl(ALPHABETSIZE,length-1-i)) {
      ret[i]=rest/power(ALPHABETSIZE,length-1-i);
      rest-=ret[i]*power(ALPHABETSIZE,length-1-i);
  //  }
  }
  return;
}

void getAssignmentFromComplementaryIndex(int index, int length, int *ret) {
  int buf[length];
  int i;

  getAssignmentFromIndex(index, length, buf);
  for (i=0; i<length; i++) {
    ret[length-1-i]=getComplementFromIndex(buf[i]);
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
    index+=getComplementNucIndex(ass[length-1-i]) *power(ALPHABETSIZE,length-i-1);
  }
  return index;
}

