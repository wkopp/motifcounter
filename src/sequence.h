#ifndef sequence_h
#define sequence_h
#include "stdio.h"

#define ALPHABETSIZE 4

typedef struct {
    char **seq;
    int  nseq;
    int  *lseq;
} Sequence;

void allocSequence(Sequence *seq, int nseq, int *lseq);
void destroySequence(Sequence *seq);
int getSequence(FILE *f,Sequence *);
int getNucIndex(char c);
char getNuc(int i);
void writeSequenceBinary(FILE *f, Sequence *s, int dhssize);
int getSequenceFromBinary(FILE *f,Sequence *s);
int getComplementFromIndex(int i);
int getIndex(int i);
void getAssignmentFromIndex(int index, int length, int *ret);
int getIndexFromReverseAssignment(char *ass, int length);
void getAssignmentFromComplementaryIndex(int index, int length, int *ret);
int getIndexFromComplementaryAssignment(char *ass, int length);
int getIndexFromAssignment(char *ass, int length);
int getIndexFromReverseComplementaryAssignment(char *ass, int length);
int skipAssignment(char *ass, int len);
int isNucleotide(char c);

#endif
