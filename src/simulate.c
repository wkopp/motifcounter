#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "forground.h"
#include "background.h"
#include "countdist.h"
#include "sequence.h"
#include "scorefunctions.h"
#include "minmaxscore.h"
#include "score1d.h"

char sampleNucleotide(double *prob) {
#ifndef IN_R
  int ip=rand();
#else
  double ip=unif_rand();
#endif
  double p;
  double cumprob=0.0;
  char n;
  int i;

  p=(double)(ip);
  #ifndef IN_R
  p/=(double)RAND_MAX;
  #endif
  for (i=0; i<ALPHABETSIZE; i++) {
    cumprob+=prob[i];
    if(p<=cumprob) {
      break;
    }
  }
  if (i==4)
    n= getNuc(i-1);
  else
    n= getNuc(i);
  return n;
}

void sampleInitialNucleotide(double *prob, char *seq, int order) {
#ifndef IN_R
  int ip=rand();
#else
  double ip=unif_rand();
#endif
  double p;
  double cumprob=0.0;
  int i, ass[order];

  p=(double)(ip);
  #ifndef IN_R
  p/=(double)RAND_MAX;
  #endif
  for (i=0; i<power(ALPHABETSIZE, order); i++) {
    cumprob+=prob[i];
    if(p<=cumprob) {
      break;
    }
  }
  getAssignmentFromIndex(i, order, ass);
  for(i=0;i<order; i++) {
    seq[i]=getNuc(ass[i]);
  }
}

void generateRandomSequence(double *station, double *trans, char *seq, int seqlen, int order) {
  int i, ind;

 sampleInitialNucleotide(station, seq, order);
 for (i=order; i<seqlen; i++) {
   ind=getIndexFromAssignment(&seq[i-order], order);
   seq[i]=sampleNucleotide(&trans[ind*ALPHABETSIZE]);
 }
}

void randomStatistics(char *seq, int seqlen, double *stat, double *trans) {
  int i;

  stat[getNucIndex(seq[0])]++;
  for (i=1; i<seqlen; i++) {
    stat[getNucIndex(seq[i])]++;
    trans[getNucIndex(seq[i-1])*ALPHABETSIZE+getNucIndex(seq[i])]++;
  }
}
void normalizeStatistics(double *stat, double *trans) {
  int i, j, sum=0;
  for (i=0; i<power(ALPHABETSIZE, 1);i++) {
  sum=0;
    for (j=0; j<power(ALPHABETSIZE, 1); j++) {
      sum+=trans[i*ALPHABETSIZE +j];
    }
    for (j=0; j<power(ALPHABETSIZE, 1); j++) {
      trans[i*ALPHABETSIZE +j]/=sum;
    }
  }
  sum=0;
  for (i=0; i<power(ALPHABETSIZE, 1);i++) sum+=stat[i];
  for (i=0; i<power(ALPHABETSIZE, 1);i++) stat[i]/=sum;
}

int countOccurances(double *station, double *trans, 
  DMatrix *pwm, DMatrix *cpwm, char *seq, int seqlen,
 int qalpha, double granularity, int order) {
 int count=0, palhits=0;
 int f=0, jump;
 int index;
 int s, i,j;
 int score[power(ALPHABETSIZE, order+1)];

 for (i=0; i< seqlen-pwm->nrow+1; i++) {
   // do not calculate the score if there are non-nucleotides in the sequence
   jump=0;
   for (j=0; j<pwm->nrow; j++) {
     if (getNucIndex(seq[i+j])<0) {
       jump=1; break;
     }
   }
   if(jump>0) continue;


   s=0;
   index=0;
   if (order>0) {
     getScoresInitialIndex(pwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s+=score[index];
   }
   for (j=order; j<pwm->nrow; j++) {

     index*=ALPHABETSIZE;
     index+=getNucIndex(seq[i+j]);

     s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index], granularity);

     index-=(index/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
   }
   if (s>=qalpha) {
     count++;
     f=1;
   }


   s=0;
   index=0;
   if (order>0) {
     getScoresInitialIndex(cpwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s+=score[index];
   }

   for (j=order; j<pwm->nrow; j++) {

     index*=ALPHABETSIZE;
     index+=getNucIndex(seq[i+j]);

     s+=getScoreIndex(cpwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index], granularity);

     index-=(index/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
   }

   if (f==0&&s>=qalpha) {
     count++;
   }
   if (f==1&&s>=qalpha) {
     palhits++;
   }
   f=0;
 }
 #ifdef DEBUG
 #ifdef IN_R
 Rprintf("%d\n", count+ palhits);
 #else
 printf("%d\n", count+ palhits);
 #endif
 #endif
 return count+palhits;
}

int countOccurancesSingleStranded(double *station, double *trans, 
  DMatrix *pwm, DMatrix *cpwm, char *seq, int seqlen,
 int qalpha, double granularity, int order) {
 int count=0, palhits=0;
 //int f=0; 
 int jump;
 int index;
 int s, i,j;
 int score[power(ALPHABETSIZE, order+1)];

 for (i=0; i< seqlen-pwm->nrow+1; i++) {
   // do not calculate the score if there are non-nucleotides in the sequence
   jump=0;
   for (j=0; j<pwm->nrow; j++) {
     if (getNucIndex(seq[i+j])<0) {
       jump=1; break;
     }
   }
   if(jump>0) continue;


   s=0;
   index=0;
   if (order>0) {
     getScoresInitialIndex(pwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s+=score[index];
   }
   for (j=order; j<pwm->nrow; j++) {

     index*=ALPHABETSIZE;
     index+=getNucIndex(seq[i+j]);

     s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index], granularity);

     index-=(index/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
   }
   if (s>=qalpha) {
     count++;
     //f=1;
   }


#ifdef NOT_NEEDED_FOR_SINGLESTRANDED_CASE
   s=0;
   index=0;
   if (order>0) {
     getScoresInitialIndex(cpwm->data, station, score, &granularity, order);
     index=getIndexFromAssignment(&seq[i], order);
     s+=score[index];
   }

   for (j=order; j<pwm->nrow; j++) {

     index*=ALPHABETSIZE;
     index+=getNucIndex(seq[i+j]);

     s+=getScoreIndex(cpwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index], granularity);

     index-=(index/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
   }

   if (f==0&&s>=qalpha) {
     count++;
   }
   if (f==1&&s>=qalpha) {
     palhits++;
   }
   f=0;
#endif
 }
 #ifdef DEBUG
 #ifdef IN_R
 Rprintf("%d\n", count+ palhits);
 #else
 printf("%d\n", count+ palhits);
 #endif
 #endif
 return count+palhits;
}

void scoreOccurances(double *station, double *trans, 
  DMatrix *pwm, char *seq, int seqlen,
 double *dist, double granularity, int smin, int order) {
 int i, j;
 int is, index, jump;
 double s;


 for (i=0; i< seqlen-pwm->nrow+1; i++) {
   s=0.0;
   index=0;
   // jump section has the purpose of avoiding
   // to process any non-nucleotide letters
   jump=0;
   for (j=0; j<pwm->nrow; j++) {
     if (getNucIndex(seq[i+j])<0) {
       jump=1; break;
      }
    }
    if(jump>0) continue;

   if (order>0)  {
     for (j=0; j<order; j++) {
       s+=log(pwm->data[getNucIndex(seq[i+j])]);
       index+=getNucIndex(seq[j+i])*power(ALPHABETSIZE, order-j-1);
     }
     s-=log(station[index]);
   }

   is=(int)roundl(s/granularity);
   for (j=order; j<pwm->nrow; j++) {

     index=index*ALPHABETSIZE + getNucIndex(seq[i+j]);

     is+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
     trans[index],granularity);

     index-=(index/power(ALPHABETSIZE,order))*power(ALPHABETSIZE,order);
   }
   //is=s;
   dist[is-smin]++;
 }
}


#ifndef IN_R
void simulateCountDistribution(char *bgfile, char *pwmfile, 
   char *output_dist, char *sdx, char * spval, char* perm, 
   char *slen,  char *mxhit, char *snos) {
  int Nhits,seqlen,maxhits,Nperm, intervalsize;
  DMatrix pwm, cpwm;
  ExtremalScore escore;
  MotifScore1d null;
  double *station, *trans;
  double dx, quantile,pvalue;
  double *dist;
  int n, nos, s;
  char *seq;
  FILE *f, *fcpd;
  int order, threshold;

#ifndef IN_R
  srand(time(0));
#else
  GetRNGstate();
#endif
  seqlen=atol(slen);
  maxhits=atol(mxhit);
  Nperm=atol(perm);
  pvalue=atof(spval);
  dx=atof(sdx);
  nos=atol(snos);

  f =fopen(bgfile,"rb");
  readBackground(f,&station, &trans, &order);
  fclose(f);
  f =fopen(pwmfile,"rb");
  readMatrix(f,&pwm);
  getComplementaryMotif(&pwm,&cpwm);
  fclose(f);

  // compute significance threshold
  initExtremalScore(&escore, dx, pwm.nrow, order);

  loadMinMaxScores(&pwm, station, trans, &escore);
  loadIntervalSize(&escore, NULL);
  intervalsize=maxScoreIntervalSize(&escore);

  initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
  null.meta.prob=&ProbBg;
  null.meta.probinit=&ProbinitBg;
  initScoreDistribution1d(&pwm,trans,&null, order);

  computeScoreDistribution1d(&pwm, trans,  station, &null, &escore, order);

  quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,pvalue));
  threshold=(int)(quantile/dx);
  #ifndef IN_R
  fprintf(stdout,"threshold=%d\n",threshold);
  #endif
  deleteExtremalScore(&escore);
  deleteScoreDistribution1d(&null, order);
  // end  of significance threshold computation

  dist=initCountDistribution(maxhits);
  #ifdef IN_R
  seq=Calloc(seqlen, char);
  #else
  seq=calloc(seqlen, sizeof(char));
  #endif

  for (n=0; n<Nperm; n++) {
    Nhits=0;
    for (s=0;s<nos; s++) {
      generateRandomSequence(station, trans, seq, seqlen, order);
      Nhits+=countOccurances(station, trans, &pwm, &cpwm, seq, seqlen, threshold, dx, order);

    }
    if (Nhits>maxhits) Nhits=maxhits;
    dist[Nhits]+=1.0/(double)Nperm;
  }

  fcpd=fopen(output_dist,"w");
  writeCountDistribution(fcpd, dist, maxhits);
  fclose(fcpd);
  deleteBackground(station, trans);

  deleteMatrix(&pwm);
  deleteMatrix(&cpwm);
  deleteCountDistribution(dist);
  #ifdef IN_R
  Free(seq);
  PutRNGstate();
  #else
  free(seq);
  #endif
}

void simulateScores(char *bgmodel, char *pwmfile, char *output, char *slen, 
  char *perm, char *gran) {
  int i, n;
  DMatrix pwm;
  double *station, *trans;
  double *estat, *etrans;
  FILE *f, *fout;
  char *seq, nuc[]="acgt";
  ExtremalScore escore;
  int mins, maxs, j;
  int seqlen, Nperm, order;
  double dx, *dist;

#ifndef IN_R
  srand(time(0));
#else
  GetRNGstate();
#endif

  seqlen=atol(slen);
  Nperm=atol(perm);

  dx=atof(gran);

  //fprintf(stdout, "Reading %s ...\n",bgmodel);
  f =fopen(bgmodel,"rb");
  if (f==0) {
    fprintf(stderr, "Error during reading %s\n",bgmodel);
    return;
  }
  readBackground(f,&station, &trans, &order);
  fclose(f);
  f=NULL;
  //fprintf(stdout, "Reading %s ...\n",pwmfile);
  f =fopen(pwmfile,"rb");
  if (f==0) {
    fprintf(stderr, "Error during reading %s\n",pwmfile);
    return;
  }
  readMatrix(f,&pwm);
  fclose(f);
    //seqlen=pwm.nrow;
  #ifdef IN_R
  estat=Calloc(power(ALPHABETSIZE,order+1), double);
  etrans=Calloc(power(ALPHABETSIZE,order+1), double);

  seq=Calloc(seqlen, char);
  #else
  estat=calloc(power(ALPHABETSIZE,order+1), sizeof(double));
  etrans=calloc(power(ALPHABETSIZE,order+1), sizeof(double));

  seq=calloc(seqlen, sizeof(char));
  #endif


  initExtremalScore(&escore, dx, pwm.nrow, order);

  loadMinMaxScores(&pwm, station, trans, &escore);
  loadIntervalSize(&escore, NULL);

  mins=getTotalScoreLowerBound(&escore);
  maxs=getTotalScoreUpperBound(&escore);
        //initScoreMetaInfo(maxs,mins,dx, &alterbf.meta);
        //initScoreMetaInfo(maxs,mins,dx, &nullbf.meta);

#ifdef IN_R
  dist=Calloc(maxs-mins+1, double);
  #else
  dist=calloc(maxs-mins+1, sizeof(double));
  #endif

  for (n=0; n<Nperm; n++) {
    generateRandomSequence(station, trans, seq, seqlen, order);
    if(order<2) randomStatistics(seq, seqlen, estat, etrans);
    scoreOccurances(station, trans, &pwm, seq, seqlen, dist, dx, mins,order);
  }
  if (order<2) {
    normalizeStatistics(estat, etrans);
    for (i=0; i<power(ALPHABETSIZE,order); i++) {
      for (j=0; j<ALPHABETSIZE; j++) {
        fprintf(stdout,"e(%c|%c)=%f\n", nuc[j],nuc[i],trans[i*ALPHABETSIZE+j]-etrans[i*ALPHABETSIZE+j]);
      }
    }
    for (i=0;i<4;i++) {
      fprintf(stdout,"e(%c)=%f\n", nuc[i],station[i]-estat[i]);
    }
  }
  #ifdef IN_R
  Free(estat);Free(etrans);
  #else
  free(estat);free(etrans);
  #endif
  for (n=0; n<maxs-mins+1; n++) {
    dist[n]/=(double)Nperm*(seqlen-pwm.nrow+1);
  }

  fout=fopen(output,"w");
  //writeCountDistribution(fcpd, dist, maxhits);
  for (i=0; i<maxs-mins+1; i++) {
    fprintf(fout, "%e\t%e\n", (double)(mins+i)*dx, dist[i]);
  }
  fclose(fout);

  #ifdef IN_R
  Free(dist);
  PutRNGstate();
  #else
  free(dist);
  #endif
}
#endif
