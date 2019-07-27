#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include "sequence.h"
#include "matrix.h"

extern int Rnucmap[], Rrevcomp[];

int getNucIndex(char c) {
    int l = toupper(c) - 'A';
    return Rnucmap[l];
}

int getComplementNucIndex(char c) {
    int l = toupper(c) - 'A';
    return Rrevcomp[Rnucmap[l]];
}

char getNuc(int i) {
    char r = 'n';
    switch (i) {
    case 0:
        r = 'a';
        break;
    case 1:
        r = 'c';
        break;
    case 2:
        r = 'g';
        break;
    case 3:
        r = 't';
        break;
    }
    return r;
}


void getAssignmentFromIndex(int index, int length, int *ret) {
    int rest = index;
    int i;

    for (i = 0; i < length; i++) {
        ret[i] = rest / power(ALPHABETSIZE, length - 1 - i);
        rest -= ret[i] * power(ALPHABETSIZE, length - 1 - i);
    }
    return;
}

int getIndexFromAssignment(const char *ass, int length) {
    int index = 0;
    int i;

    for (i = 0; i < length; i++) {
        index += getNucIndex(ass[i]) * power(ALPHABETSIZE, length - 1 - i);
    }
    return index;
}

int getIndexFromReverseAssignment(const char *ass, int length) {
    int index = 0;
    int i;

    for (i = 0; i < length; i++) {
        index += getNucIndex(ass[length - 1 - i]) * power(ALPHABETSIZE,
                 length - 1 - i);
    }
    return index;
}

int getIndexFromComplementaryAssignment(const char *ass, int length) {
    int index = 0;
    int i;

    for (i = 0; i < length; i++) {
        index += getComplementNucIndex(ass[i]) * power(ALPHABETSIZE, length - 1 - i);
    }
    return index;
}

int getIndexFromReverseComplementaryAssignment(const char *ass, int length) {
    int index = 0;
    int i;

    for (i = 0; i < length; i++) {
        index += getComplementNucIndex(ass[length - 1 - i]) *
                 power(ALPHABETSIZE, length - i - 1);
    }
    return index;
}

int getSequenceLength(const char *seq, int slen) {
    return slen;
}

int hasN(const char *seq, int slen) {
  int j;
  
  for (j = 0; j < slen; j++) {
    // if the sequence contains any N's, do not process the scores
    if (getNucIndex(seq[j]) < 0) {
      return 1;
    }
  }
  return 0;
}
