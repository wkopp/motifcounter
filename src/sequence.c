#include <stdio.h>
#include <R.h>
#include "sequence.h"
#include "matrix.h"

int getNucIndex(char c) {
    int r = -1;
    switch (c) {
    case 'a':
    case 'A':
        r = 0;
        break;
    case 'c':
    case 'C':
        r = 1;
        break;
    case 'g':
    case 'G':
        r = 2;
        break;
    case 't':
    case 'T':
        r = 3;
        break;
    }
    return r;
}

int getComplementNucIndex(char c) {
    int r = -1;
    switch (c) {
    case 'a':
    case 'A':
        r = 3;
        break;
    case 'c':
    case 'C':
        r = 2;
        break;
    case 'g':
    case 'G':
        r = 1;
        break;
    case 't':
    case 'T':
        r = 0;
        break;
    case 'n':
    case 'N':
        r = -1;
        break;
    }
    return r;
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
    int j;

    for (j = 0; j < slen; j++) {
        // if the sequence contains any N's, do not process the scores
        if (getNucIndex(seq[j]) < 0) {
            return -1;
        }
    }
    return slen;
}
