
#ifndef background_h
#define background_h


void getNucleotideFrequencyFromSequence(char *seq, int slen,
        double *counts, int order);
int getForwardTransition(double *di, double *forwardtrans, int order);

void readBackground (FILE *f, double **station, double **trans, int *order);
void printBackground(double *stat, double *trans, int);
void deleteBackground(double * station, double *trans);
int getStationaryDistribution(double *trans, double *station, int order);
int makebg(char *infasta, char *outmodel, char *outseq, char *order);
void writeBackground(FILE *f, double * station, double * trans, int order);
#endif
