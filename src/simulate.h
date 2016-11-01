#ifndef simulate_h
#define simulate_h

void generateRandomSequence(double *station, double *trans, char *seq, int seqlen, int order);
void randomStatistics(char *seq, int seqlen, double *stat, double *trans);
int countOccurances(double *station, double *trans, DMatrix *pwm, char *seq, int seqlen,
 int qalpha, double granularity, int order);
void scoreOccurances(double *station, double *trans, 
  DMatrix *pwm, char *seq, int seqlen,
 double *dist, double granularity, int smin, int order);
void simulateCountDistribution(char *bgfile, char *pwmfile, 
   char *output_dist, char *sdx, char * spval, char* perm, 
   char *slen,  char *mxhit, char *snos);
//void simulateCountDistribution(char *bgfile, char *pwmfile, char *output_quantile, 
//   char *output_dist, char *gr, double alphaquantile, int Nperm, int seqlen);
void simulateScores(char *bgmodel, char *pwmfile, char *output, char *seqlen, 
  char *perm, char *gran);
void randomStatistics(char *seq, int seqlen, double *stat, double *trans);
void normalizeStatistics(double *stat, double *trans);

#endif
