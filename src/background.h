
#ifndef background_h
#define background_h

#include "matrix.h"

//void getBackgroundFromSequence(FILE *f, DMatrix *m);
//int getBackgroundFromSequence(FILE *f,DMatrix *mono, DMatrix *di);
//int getOrder1BackgroundFromSequence(FILE *f, DMatrix *di);
int getNucleotideFrequencyFromSequence(FILE *f, double *di, int order);
int getForwardTransition(double *di, double *forwardtrans, int order);
//int getReverseTransition(DMatrix *di, DMatrix *forwardtrans);
//int getReverseTransition(DMatrix *di, DMatrix *stationary, DMatrix *reversetrans);
void readBackground (FILE *f, double **station, double **trans, int *order);
void deleteBackground(double * station, double *trans);
int getStationaryDistribution(double *trans, double *station, int order);
int makebg(char *infasta, char *outmodel, char *outseq, char *order);
void writeBackground(FILE *f, double * station, double * trans, int order);
#endif
