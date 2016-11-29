
#ifndef sharedtypes_h
#define sharedtypes_h

#include <stdio.h>
#include <stdlib.h>


typedef struct {
    int nrow,ncol;
    double *data;
} DMatrix;

typedef struct {
    int nrow,ncol;
    int *data;
} IMatrix;

int power(int base, int exp);
#endif
