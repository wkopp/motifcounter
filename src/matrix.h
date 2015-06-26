
#ifndef sharedtypes_h
#define sharedtypes_h

#include <stdio.h>
#include <stdlib.h>

//#define WSIZE 550
//#define BSIZE 100*550
//#define MAXMOTIFWIDTH 151
//#define MAXBGWIDTH 100000+MAXMOTIFWIDTH
#define DHSSIZE 150

typedef struct {
  int nrow,ncol;
  double *data;
} DMatrix;

typedef struct {
  int nrow,ncol;
  int *data;
} IMatrix;

void deleteMatrix(DMatrix *m);
int writeMatrix(FILE *f, DMatrix *m);
int readMatrix(FILE *f,DMatrix *m);
int power(int base, int exp);
#endif
