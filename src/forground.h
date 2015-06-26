
#include "matrix.h"
int getMotifWidth(FILE *f);
void getMotif(FILE *f, DMatrix *m, double);
void getComplementaryMotif(DMatrix *pwm, DMatrix *cpwm);
int makemotif(char *fmotif, char *binmotif, char*pseudo, char *format);
void printMotif(DMatrix *pwm);
int getTableMotifWidth(FILE *f);
int getTransfacMotifWidth(FILE *f);
int getJasparMotifWidth(FILE *f);
void getTableMotif(FILE *f, DMatrix *m, double pseudocount);
void getTransfacMotif(FILE *f, DMatrix *m, double pseudocount);
void getJasparMotif(FILE *f, DMatrix *m, double pseudocount);
