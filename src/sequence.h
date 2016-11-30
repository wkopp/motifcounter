#ifndef sequence_h
#define sequence_h
#include "stdio.h"

#define ALPHABETSIZE 4

int getNucIndex(char c);
char getNuc(int i);
void getAssignmentFromIndex(int index, int length, int *ret);
int getIndexFromReverseAssignment(char *ass, int length);
int getIndexFromComplementaryAssignment(char *ass, int length);
int getIndexFromAssignment(char *ass, int length);
int getIndexFromReverseComplementaryAssignment(char *ass, int length);

#endif
