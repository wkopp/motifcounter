#ifndef sequence_h
#define sequence_h
#include "stdio.h"

#define ALPHABETSIZE 4

int getNucIndex(char c);
char getNuc(int i);
void getAssignmentFromIndex(int index, int length, int *ret);
int getIndexFromReverseAssignment(const char *ass, int length);
int getIndexFromComplementaryAssignment(const char *ass, int length);
int getIndexFromAssignment(const char *ass, int length);
int getIndexFromReverseComplementaryAssignment(const char *ass, int length);
int getSequenceLength(const char *seq, int slen);
int hasN(const char *seq, int slen);

#endif
