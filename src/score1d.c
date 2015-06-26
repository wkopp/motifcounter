#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "scorefunctions.h"
#include "score1d.h"
#include "forground.h"
#include "compoundpoisson.h"
#include "countdist.h"
#include "overlap.h"
#include "background.h"
//#include "inputoutput.h"

#undef PARALLEL_WK
#define PARALLEL_WK

double getQuantileWithIndex1d(MotifScore1d *s, int qi) {
  //fprintf(stdout,"xmin=%d\n",s->meta.xmin);
  return (double)(s->meta.xmin + qi)*s->meta.dx;
}

int getQuantileIndex1d(Score1d *s, double pvalue) {
  int i, x=1;
  double p=0;
  for (i=s->end; i>=0; i--) {
    p+=s->y[i];
    if (pvalue<p) { break;}
    if(s->y[i]==0.0) {
      x++;
    } else {
      x=1;
    }
  //  fprintf(stdout, "p=%f\n",p);
  }
//  fprintf(stdout,"i=%d, x=%d p=%f, alpha=%f\n",i, x,p,pvalue);
  return i+x;
}

double getProbWithIndex1d(MotifScore1d *s, int iquantile) {
  int i;
  double sum=0.0;

  for (i=iquantile-s->meta.xmin; i<=s->totalScore.end; i++) {
    sum+=s->totalScore.y[i];
  }
 // fprintf(stdout,"iquantile=%d, xmin=%d,end=%d, sum=%f\n",
 // iquantile,s->meta.xmin,s->totalScore.end,sum);
  return sum;
}

double getProb1d(MotifScore1d *s, double quantile) {
  double r=(quantile-s->meta.xmin)/(s->meta.dx);
  int iq=(int)roundl(r);

  return getProbWithIndex1d(s,iq);
}

void storeScoreDist1d (FILE *f, MotifScore1d *s, int withhead) {
  int i;
  if (withhead==1) {
  for (i=0; i<=s->meta.xmax-s->meta.xmin; i++) {
    fprintf(f, "%e ", (double)(s->meta.xmin+i)*s->meta.dx);
  }
  fprintf(f, "\n");
  }
  for (i=0; i<=s->meta.xmax-s->meta.xmin; i++) {
  //for (i=0; i<=s->totalScore.end-s->totalScore.start; i++) {
    fprintf(f, "%e ", s->totalScore.y[i]);
  }
  fprintf(f, "\n");
}


void initScore1d(Score1d *s, int l) {
#ifdef IN_R
  s->y=Calloc(l, double);
  #else
  s->y=calloc(l, sizeof(double));
  #endif
  s->end=0;
  s->merged=0;
  s->start=l;
}

int initScoreDistribution1d (DMatrix *theta, double *bg1, MotifScore1d *result, int order) {
  int i;

  initScore1d(&result->totalScore, result->meta.xmax-result->meta.xmin+1);

  result->mlen=theta->nrow;
  #ifdef IN_R
  result->ScoreBuffer1=Calloc(power(ALPHABETSIZE, order)*theta->nrow, Score1d);
  result->tmpScore=Calloc(power(ALPHABETSIZE, order+1), Score1d);
  #else
  result->ScoreBuffer1=calloc(power(ALPHABETSIZE, order)*theta->nrow, sizeof(Score1d));
  result->tmpScore=calloc(power(ALPHABETSIZE, order+1), sizeof(Score1d));
  #endif

  for (i=0; i < power(ALPHABETSIZE, order)*theta->nrow; i++) {
    initScore1d(&result->ScoreBuffer1[i], result->meta.length);
  }

  for (i=0; i < power(ALPHABETSIZE, order+1); i++) {
    initScore1d(&result->tmpScore[i], result->meta.length);
  }
  return 0;
}

int deleteScoreDistribution1d(MotifScore1d *m, int order) {
  int j;

#ifdef IN_R
  for (j=0; j < power(ALPHABETSIZE, order)*m->mlen; j++) {
    Free(m->ScoreBuffer1[j].y);
  }
  for (j=0; j < power(ALPHABETSIZE, order+1); j++) {
    Free(m->tmpScore[j].y);
  }

  Free(m->ScoreBuffer1);
  Free(m->tmpScore);

  Free(m->totalScore.y);
  #else
  for (j=0; j < power(ALPHABETSIZE, order)*m->mlen; j++) {
    free(m->ScoreBuffer1[j].y);
  }
  for (j=0; j < power(ALPHABETSIZE, order+1); j++) {
    free(m->tmpScore[j].y);
  }

  free(m->ScoreBuffer1);
  free(m->tmpScore);

  free(m->totalScore.y);
  #endif
  return 0;
}

void resetScore1d(Score1d *score, ScoreMetaInfo *meta) {
  memset(score->y,0,(meta->length)*sizeof(double));
  score->start=meta->length;
  score->end=0;
}

void addScore1d(Score1d *a, Score1d *b, ScoreMetaInfo *meta) {
  int i;
  if (b->start > b->end) return;
  a->start=(a->start<b->start) ? a->start : b->start;
  a->end=(a->end>b->end) ? a->end : b->end;
  #ifdef PARALLEL_WK
  #pragma omp parallel for default(none) shared(a,b) private(i)
  #endif
  for (i=b->start; i<=b->end;i++) {
    a->y[i]+=b->y[i];
  }
}

double getProbability1d(Score1d *a, ScoreMetaInfo *meta) {
  int i;
  double sum=0;
  for (i=a->start; i<=a->end;i++){
    sum+=a->y[i];
  }
  return sum;
}

#ifdef WKOPP
void normalizeScore1d(Score1d *a, ScoreMetaInfo *meta) {
  int i;
  double sum=0;
  for (i=0; i<meta->length;i++){
    sum+=a->y[i];
  }
  for (i=0; i<meta->length;i++){
    a->y[i]/=sum;
  }
}
#endif

static void ShiftMultiplyScoreIndex1d(Score1d *dest, Score1d *src, 
  int *ds, double p, int lprev, int uprev, int lcur, int ucur) {
  int i;
  int restinterval;

  
  if (src->start > src->end) return;
  if (p == 0.0) return;

  // lprev and lcur are used to reconstruct the absolute index, since we are only
  // using relative indexes with src->start and src->end
  //
  // compute putative dest->start
  // if dest->start < 0, set start to 0. probability mass is dropped out then, because
  // it cannot reach the threshold anymore.
  // if dest->end > intervalsize at this position, set dest->end to lcur + intervalsize.
  // probability mass beyond that cannot deceed the threshold.
  // use memmove to copy the values between dest->start and dest->end.
  // add up values that cannot deceed the threshold.
  //
  // finally multiply all values in dest->y by p.
  //fprintf(stderr, "s.s=%d, s.e=%d, lp=%d, up=%d\n lc=%d,uc=%d\n",
 // src->start, src->end, lprev,uprev,lcur,ucur);

  dest->start = src->start-lcur+lprev+*ds;
  if (dest->start<0) {
    src->start-=dest->start;
    dest->start=0;
  }

  dest->end = src->end-lcur+lprev+*ds;
  if (dest->end > ucur -lcur) {
    restinterval=dest->end- ucur+lcur;
    dest->end=ucur -lcur;
  } else {
    restinterval=0;
  }

  //fprintf(stderr, "s.s=%d, s.e=%d, lp=%d, up=%d\n lc=%d,uc=%d\n",
  //src->start, src->end, lprev,uprev,lcur,ucur);
  //fprintf(stderr, "e.s=%d, e.e=%d, resti=%d\n", dest->start,dest->end,restinterval);

  if (dest->start>dest->end) { return; }

  //dest->end = (dest->end>ucur) ? ucur : dest->end;
  memmove(&dest->y[dest->start],&src->y[src->start], (dest->end-dest->start+1)*sizeof(double));

  #ifdef PARALLEL_WK
  #pragma omp parallel for default(none) shared(p,dest) private(i)
  #endif
  for (i=dest->start; i<=dest->end;i++) {
    dest->y[i]*=p;
  }
  for (i=0; i<restinterval;i++) {
    dest->y[dest->end]+=src->y[src->start+dest->end-dest->start+1+i]*p;
  }
}

int computeTotalScoreDistribution1d(MotifScore1d *mscore, ExtremalScore *tm, int order) {

  int i=0, k;
  int lmin=getTotalScoreLowerBound(tm);
  int *lbound=getLastScoreLowerBound(tm);
  //int L[]={343,325,319,233};
  //double P[4];
  Score1d *a,*b;
  a=&mscore->totalScore;
  for (i=0; i<power(ALPHABETSIZE,order); i++) {
   // P[i]=0.0;
    b=&mscore->ScoreBuffer1[power(ALPHABETSIZE,order)*(mscore->mlen-1) + i];

    if (b->start>b->end) return 0;
    a->start=(a->start<b->start+lbound[i]-lmin) ? a->start : b->start+lbound[i]-lmin;
    a->end=(a->end>b->end+lbound[i]-lmin) ? a->end : b->end+lbound[i]-lmin;
  #ifdef PARALLEL_WK
    #pragma omp parallel for default(none) shared(a,b,lbound,lmin, i) private(k)
    #endif
    for (k=b->start; k<=b->end;k++) {
      a->y[k+lbound[i]-lmin]+=b->y[k];
    }
  }
  return 0;
}

#ifdef WKO
int getIntervalSize(DMatrix *pwm1, DMatrix *pwm2, double *station, double *trans, int order,ExtremalScore *es1,ExtremalScore *es2) {
  int size1, size2;
    // 1D score
    loadMinMaxScores(pwm1, station, trans, es1);
    loadMinMaxScores(pwm2, station, trans, es2);
    loadIntervalSize(es1,NULL);
    loadIntervalSize(es2,NULL);
    size1=maxScoreIntervalSize(es1);
    size2=maxScoreIntervalSize(es2);

    return (size1>size2) ? size1 : size2;
}

void computeScoreDist1d(DMatrix *pwm1, DMatrix *pwm2, double *station, double *trans, double *dx
  int order, MotifScore1d *mscore1, MotifScore1d *mscore2) {
  ExtremalScore es1, es2;
  //MotifScore1d init1d1, init1d2;
  initExtremalScore(&es1, *dx, pwm1->nrow, order);
  initExtremalScore(&es2, *dx, pwm1->nrow, order);

  intervalsize=getIntervalSize(pwm1,pwm2,station, trans, order, &es1,&es2);

  initScoreMetaInfo(getTotalScoreLowerBound(&es1),getTotalScoreUpperBound(&es2),
    intervalsize, *dx, mscore1->meta);
  initScoreMetaInfo(getTotalScoreLowerBound(&es1),getTotalScoreUpperBound(&es2),
    intervalsize, *dx, mscore2->meta);

    initScoreDistribution1d(pwm1,trans,mscore1, order);
    initScoreDistribution1d(pwm2,trans,mscore2, order);

    computeScoreDistribution1d(pwm1, trans,  station, mscore1, &es1, order);
    computeScoreDistribution1d(pwm2, trans,  station, mscore2, &es2, order);

    deleteExtremalScore(&es1);
    deleteExtremalScore(&es2);
}

int computeScoreThresholdFromPvalue(DMatrix *pwm1, DMatrix *pwm2, double *station, double *trans,
  double *pvalue, double *dx, int order) {
  MotifScore1d mscore1, mscore2;

    computeScoreDist1d(pwm1,pwm2,station, trans,dx,order,&mscore1,&mscore2);

    quantile=getQuantileWithIndex1d(&mscore1,getQuantileIndex1d(&mscore1.totalScore,*pvalue));
    threshold=(int)(quantile/(*dx));
    #ifdef WK
    fprintf(stdout, "Probability mass-null: %f\n",
    getProbability1d(&init1d1.totalScore, &init1d1.meta));
    fprintf(stdout, "pval=%f, th=%f th=%d pvalf=%f pvalr=%f\n", *pvalue, quantile,threshold,
    getProbWithIndex1d(&init1d1, threshold),
    getProbWithIndex1d(&init1d2, threshold));
    #endif

    deleteScoreDistribution1d(&mscore1, order);
    deleteScoreDistribution1d(&mscore2, order);
    //loadIntervalSize(&uescore1, &threshold);
    //loadIntervalSize(&uescore2, &threshold);
    //cutScoreRangeWithThreshold(&init1d1, &uescore1, order);
    //cutScoreRangeWithThreshold(&init1d2, &uescore2, order);

  return threshold;
}
#endif
#define DEBUG
#undef DEBUG
void cutScoreRangeWithThreshold(MotifScore1d *mscore, ExtremalScore *tm, int order) {
  int m, i, l;
  int start, end, rest;
//  double p1=0.0, p2=0.0;
  if (order==0) {
    m=0;
  } else {
    m=order-1;
  }

  for (; m<tm->len;m++) {
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
      if(getScoreUpperBoundUnconstrainted(tm,m,i)<getScoreLowerBound(tm, m, i)) {
        memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y, 0,
          mscore->meta.length*sizeof(double));
        //mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[0]=0;
        mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].start=1;
        mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].end=0;
        continue;
      }
      start=getScoreLowerBound(tm, m, i)-getScoreLowerBoundUnconstrainted(tm,m,i);
      end=getScoreUpperBound(tm, m, i)-getScoreLowerBoundUnconstrainted(tm,m,i);
      rest=getScoreUpperBoundUnconstrainted(tm,m,i) - getScoreUpperBound(tm, m, i);
      if (start<0) {
        start=0; end=0; rest=0;
      }
      //fprintf(stdout, "s=%d e=%d r=%d ",start,end,rest);

      mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].start=0;
      mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].end=
       getScoreUpperBound(tm,m,i)-getScoreLowerBound(tm, m, i);

      if (start>0) {
        if (mscore->meta.length>start) {
          memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y, 0,
            start*sizeof(double));
        } else {
          memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y, 0,
            mscore->meta.length*sizeof(double));
        }
      }
      for (l=0; l<=end-start; l++) {
        mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l]=
          mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l+start];
        if (start>0) {
          mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l+start]=0;
        }
      }
      for (l=0; l<rest; l++) {
        mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[end-start]+=
          mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l+end+1];
        mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l+end+1]=0;
      }
      #ifndef IN_R
      #ifdef DEBUG
      fprintf(stdout, "m=%d, i=%d, [%d,%d]\n",m, i,
        mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].start,
        mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].end);
        #endif
        #endif
    }
  }
  resetScore1d(&mscore->totalScore, &mscore->meta);
  for (i=0; i<power(ALPHABETSIZE, order); i++) {
    //fprintf(stdout,"p(%d)=%f\n",i,
    //  mscore->ScoreBuffer1[(tm->len-1)*power(ALPHABETSIZE, order) +i].y[0]);
    mscore->totalScore.y[0]+=mscore->ScoreBuffer1[(tm->len-1)*power(ALPHABETSIZE, order) +i].y[0];
  }
  mscore->totalScore.start=0;
  mscore->totalScore.end=0;
}

int computeScoreDistribution1d(DMatrix *pwm, double *trans, 
  double *station, MotifScore1d *mscore, 
  ExtremalScore *tm, int order) {
  int i,m,j, ji, k, K, corder;
  int score[power(ALPHABETSIZE, order+1)];

  if (order > pwm->nrow) {
    #ifdef IN_R
    error("Background order cannot be longer than the motif.\n");
    #else
    fprintf(stderr, "Background order cannot be longer than the motif.\n");
    #endif
    return 1;
  }
  corder=((order==0) ? (order+1) : order);
  getScoresInitialIndex(pwm->data,station, score, &mscore->meta.dx, order);


  for (i=0; i<power(ALPHABETSIZE, order); i++) {
    resetScore1d(&mscore->ScoreBuffer1[i], &mscore->meta);
    if (order==0) {
      k=0;
      K=ALPHABETSIZE;
    } else {
      k=i;
      K=i+1;
    }
    mscore->ScoreBuffer1[(corder-1)*power(ALPHABETSIZE, order) +i].start=0;
    mscore->ScoreBuffer1[(corder-1)*power(ALPHABETSIZE, order) +i].end=
        getScoreUpperBound(tm,corder-1,i)-getScoreLowerBound(tm, corder-1, i);

    for (;k<K;k++) {
      if (score[k]<getScoreLowerBound(tm, corder-1, i)) { continue; }

      mscore->ScoreBuffer1[(corder-1)*power(ALPHABETSIZE, order) +i].y[
      score[k]-getScoreLowerBound(tm,corder-1,i)]+=
      mscore->meta.probinit(station[k],pwm->data, k, corder);
    }
  }

  for (m=corder; m<pwm->nrow;m++) {
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
      getScoresIndex(&pwm->data[m*ALPHABETSIZE],&trans[i*ALPHABETSIZE], 
                score,&mscore->meta.dx);

      for (j=0; j<ALPHABETSIZE; j++) {

        resetScore1d(&mscore->tmpScore[i*ALPHABETSIZE+j], &mscore->meta);
        if (order>0) {
          ji=i*ALPHABETSIZE + j;
          ji-=(ji/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
        } else {
          ji=0;
        }

        ShiftMultiplyScoreIndex1d(&mscore->tmpScore[i*ALPHABETSIZE+j], 
             &mscore->ScoreBuffer1[(m-1)*power(ALPHABETSIZE, order)+i], 
             &score[j],
             mscore->meta.prob(trans[i*ALPHABETSIZE+j], pwm->data[m*ALPHABETSIZE +j]),
             getScoreLowerBound(tm, m-1, i),
             getScoreUpperBound(tm, m-1, i),
             getScoreLowerBound(tm, m, ji),
             getScoreUpperBound(tm, m, ji));

       // fprintf(stderr, "m=%d, i=%d, j=%d, ji=%d, lb=%d\n",m, i,j,ji,getScoreLowerBound(tm, m, ji));
      }
    }
    for (i=0; i<power(ALPHABETSIZE, order); i++) {
      //resetScore1d(&mscore->ScoreBuffer1[i], &mscore->meta);
      for (j=0; j<ALPHABETSIZE; j++) {
        addScore1d(&mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i], 
        &mscore->tmpScore[i +j*power(ALPHABETSIZE, order)], &mscore->meta);
      }
      //fprintf(stdout, "s=%d, e=%d\n",mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].start,
      // mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].end);
    }
  }
  //resetScore1d(&mscore->totalScore, &mscore->meta);
  computeTotalScoreDistribution1d(mscore, tm, order);

  return 0;
}

int computeMarginalScoreDistribution1dBruteForce(DMatrix *pwm, double *trans, 
  double *station, MotifScore1d *mscore, int xmin, int order) {
  int i, n, x, l1;
  int si, cletter, prefix, s;
  int score[power(ALPHABETSIZE, order+1)];
  double p;
  int corder=order;
  if (corder==0) corder++;

  if (order > pwm->nrow) {
    #ifdef IN_R
    error("Background order cannot be longer than the motif.\n");
    #else
    fprintf(stderr, "Background order cannot be longer than the motif.\n");
    #endif
    return 1;
  }
  getScoresInitialIndex(pwm->data,station, score, &mscore->meta.dx, order);

  for (i=0; i<power(ALPHABETSIZE,pwm->nrow); i++) {

    prefix=i/power(ALPHABETSIZE, pwm->nrow-corder);

    si=score[prefix];
    if (order>0) {
      p=mscore->meta.probinit(station[prefix], pwm->data,prefix,order);
    } else {
      p=mscore->meta.probinit(station[prefix], pwm->data,prefix,order+1);
    }

    l1=i-prefix*power(ALPHABETSIZE, pwm->nrow-corder);

    prefix=(order==0) ? 0 : prefix;
    for (n=corder; n< pwm->nrow;n++) {
      cletter=l1/power(ALPHABETSIZE, pwm->nrow-n-1);

      l1-=cletter*power(ALPHABETSIZE, pwm->nrow-n-1);
      s=getScoreIndex(pwm->data[n*ALPHABETSIZE+cletter],
        trans[prefix*ALPHABETSIZE+ cletter], mscore->meta.dx);
      p*=mscore->meta.prob(trans[prefix*ALPHABETSIZE + cletter], 
      pwm->data[n*ALPHABETSIZE + cletter]);

      si+=s;

      x=prefix*ALPHABETSIZE+cletter;
      x-=(x/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);
      prefix=x;
    }
    mscore->totalScore.y[si-xmin]+=p;
  }

  mscore->totalScore.start=0;
  mscore->totalScore.end=mscore->meta.xmax-mscore->meta.xmin;

  return 0;
}

#ifndef IN_R
void scoredist(char *bgfile, char *pwmfile, char *output_quantile, 
   char *output_dist, char *pv, char *th,char *gr) {
  int i, j,Nmotif, m;
  DMatrix *pwm, *cpwm, *motif;
  double *station, *trans;
  MotifScore1d null, alter;
  ExtremalScore escore;
  double p, quantile, dx, pcomp;
  int intervalsize;
  int order, threshold;
  FILE *f;

   if (gr==NULL) {
     fprintf(stderr, "-gran must be set\n");
     return;
  }
   dx=atof(gr);
   if (pv!=NULL) {
     p=atof(pv);
   }
   if (th!=NULL) {
     threshold=(int)atof(th)/dx;
   }
    f =fopen(bgfile,"rb");
    readBackground(f,&station, &trans, &order);
    fclose(f);
    #ifdef IN_R
    pwm=Calloc(5000,DMatrix);
    cpwm=Calloc(5000,DMatrix);
    #else
    pwm=calloc(5000,sizeof(DMatrix));
    cpwm=calloc(5000,sizeof(DMatrix));
    #endif
    f =fopen(pwmfile,"rb");
    for (i=0;i<5000;i++) {
      if(readMatrix(f,&pwm[i])!=0) {break;}
      getComplementaryMotif(&pwm[i],&cpwm[i]);
    }
    Nmotif=i;
    fclose(f);


      for (j=0; j<Nmotif; j++) {

        motif=&pwm[j];
        initExtremalScore(&escore, dx, motif->nrow, order);

        loadMinMaxScores(motif, station, trans, &escore);
        if (th!=NULL) {
          loadIntervalSize(&escore, &threshold);
        } else {
          loadIntervalSize(&escore, NULL);
        }
        intervalsize=maxScoreIntervalSize(&escore);

        initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
        initScoreMetaInfo(getTotalScoreLowerBound(&escore),
           getTotalScoreUpperBound(&escore),intervalsize,dx, &alter.meta);
        null.meta.prob=&ProbBg;
        alter.meta.prob=&ProbPWM;
        null.meta.probinit=&ProbinitBg;
        alter.meta.probinit=&ProbinitPWM;
        initScoreDistribution1d(motif,trans,&null, order);
        initScoreDistribution1d(motif,trans,&alter, order);

        computeScoreDistribution1d(motif, trans,  station, 
          &null, &escore, order);

        computeScoreDistribution1d(motif, trans,  station, 
          &alter, &escore, order);

        if (pv!=NULL) {
          quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
        } else {
          quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
        }
        threshold=(int)(quantile/dx);
        //threshold=203;
        pcomp=getProbWithIndex1d(&null,threshold);
        fprintf(stdout, "%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
        getQuantileIndex1d(&null.totalScore,p),quantile);
        if (output_dist!=NULL) {
          f=fopen(output_dist,"w");
          storeScoreDist1d(f,&null, 1);
          storeScoreDist1d(f,&alter, 0);
        }

        loadIntervalSize(&escore, &threshold);

        cutScoreRangeWithThreshold(&null, &escore, order);
        fprintf(stdout, "P(hit)= %f\n", getProbability1d(&null.totalScore,
          &null.meta));

        for (m=0;m<motif->nrow;m++) {
          for (i=0; i<power(ALPHABETSIZE, order); i++) {
            if (null.ScoreBuffer1[m*power(ALPHABETSIZE, order)+i].start>
            null.ScoreBuffer1[m*power(ALPHABETSIZE, order)+i].end) continue;
          }
        }
        if (output_dist!=NULL) {
          fclose(f);
        }
        deleteExtremalScore(&escore);
        deleteScoreDistribution1d(&null, order);
        deleteScoreDistribution1d(&alter, order);
      }
      deleteBackground(station, trans);

    for (i=0;i<Nmotif;i++) {
      deleteMatrix(&pwm[i]);
      deleteMatrix(&cpwm[i]);
    }
    #ifdef IN_R
    Free(pwm);
    Free(cpwm);
    #else
    free(pwm);
    free(cpwm);
    #endif
}

void scoredist_bruteforce(char *bgfile, char *pwmfile, char *output_quantile, 
   char *output_dist, char *pv, char *th, char *gr) {
  int i, j,Nmotif;
  DMatrix *pwm, *cpwm, *motif;
  double *station, *trans;
  MotifScore1d null, alter;
  ExtremalScore escore;
  double p, quantile, dx, pcomp;
  int intervalsize;
  int order, threshold;
  FILE *f;

   dx=atof(gr);
   if (pv!=NULL) {
     p=atof(pv);
   }
   if (th!=NULL) {
     threshold=(int)atof(th)/dx;
   }
    f =fopen(bgfile,"rb");
    readBackground(f,&station, &trans, &order);
    fclose(f);
    #ifdef IN_R
    pwm=Calloc(5000,DMatrix);
    cpwm=Calloc(5000,DMatrix);
    #else
    pwm=calloc(5000,sizeof(DMatrix));
    cpwm=calloc(5000,sizeof(DMatrix));
    #endif
    f =fopen(pwmfile,"rb");
    for (i=0;i<5000;i++) {
      if(readMatrix(f,&pwm[i])!=0) {break;}
      getComplementaryMotif(&pwm[i],&cpwm[i]);
    }
    Nmotif=i;
    fclose(f);

    for (j=0; j<Nmotif; j++) {

        motif=&pwm[j];
        initExtremalScore(&escore, dx, motif->nrow, order);

        loadMinMaxScores(motif, station, trans, &escore);
        if (th!=NULL) {
          loadIntervalSize(&escore, NULL);
        } else {
          loadIntervalSize(&escore, NULL);
        }
        intervalsize=maxScoreIntervalSize(&escore);

        initScoreMetaInfo(getTotalScoreLowerBound(&escore),
        getTotalScoreUpperBound(&escore),intervalsize,dx, &null.meta);
        initScoreMetaInfo(getTotalScoreLowerBound(&escore),
        getTotalScoreUpperBound(&escore),intervalsize,dx, &alter.meta);
        null.meta.prob=&ProbBg;
        alter.meta.prob=&ProbPWM;
        null.meta.probinit=&ProbinitBg;
        alter.meta.probinit=&ProbinitPWM;
        initScoreDistribution1d(motif,trans,&null, order);
        initScoreDistribution1d(motif,trans,&alter, order);

        computeMarginalScoreDistribution1dBruteForce(motif, trans,
          station, &null, null.meta.xmin, order);
   //     fprintf(stdout, "Probability mass-null: %f\n", getProbability1d(&null.totalScore, &null.meta));

        computeMarginalScoreDistribution1dBruteForce(motif, trans,
          station, &alter, null.meta.xmin, order);

  //      fprintf(stdout, "Probability mass-alt: %f\n", getProbability1d(&alter.totalScore, &alter.meta));
        if (pv!=NULL) {
          quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
        } else {
          quantile=getQuantileWithIndex1d(&null,getQuantileIndex1d(&null.totalScore,p));
        }
        threshold=(int)(quantile/dx);
        //threshold=203;
        pcomp=getProbWithIndex1d(&null,threshold);
        fprintf(stdout, "%1.3e-siglevel: pcomp=%1.3e, iq=%d, q=%2.2f\n",p, pcomp,
        getQuantileIndex1d(&null.totalScore,p),quantile);
        if (output_dist!=NULL) {
          f=fopen(output_dist,"w");
          storeScoreDist1d(f,&null, 1);
          storeScoreDist1d(f,&alter, 0);
          fclose(f);
        }

        threshold=(int)(quantile/dx);
        loadIntervalSize(&escore, &threshold);

        deleteExtremalScore(&escore);
        deleteScoreDistribution1d(&null, order);
        deleteScoreDistribution1d(&alter, order);
    }
    deleteBackground(station, trans);

    for (i=0;i<Nmotif;i++) {
      deleteMatrix(&pwm[i]);
      deleteMatrix(&cpwm[i]);
    }
    #ifdef IN_R
    Free(pwm);
    Free(cpwm);
    #else
    free(pwm);
    free(cpwm);
    #endif
}

void printScoreMatrix(DMatrix *pwm, double *station, double *trans, double gran, int order) {
  int i, j, m;
  int score[power(ALPHABETSIZE,order)];

  getScoresInitialIndex(pwm->data,station, score, &gran, order);
  for (i=0; i<power(ALPHABETSIZE,order); i++) {
    fprintf(stdout, "%d ", score[i]);
  }
  fprintf(stdout,"\n");
  for (m=order; m<pwm->nrow; m++) {
  fprintf(stdout,"pos:%d\n",m);
  for (i=0; i<power(ALPHABETSIZE,order); i++) {
    getScoresIndex(&pwm->data[m*ALPHABETSIZE],&trans[i*ALPHABETSIZE], score, &gran);
    for (j=0; j<ALPHABETSIZE; j++) {
      fprintf(stdout, "%d ", score[j]);
    }
    fprintf(stdout,"\n");
  }
  }
}

void printscore(char *bgfile, char *pwmfile, char *gr) {
  int i, j,Nmotif;
  DMatrix *pwm, *cpwm;
  double *station, *trans;
  double dx;
  int order;
  FILE *f;

   dx=atof(gr);
    f =fopen(bgfile,"rb");
      readBackground(f,&station, &trans, &order);
      //if(readMatrix(f,&trans[i])!=0) {break;}
    fclose(f);
    #ifdef IN_R
    pwm=Calloc(5000,DMatrix);
    cpwm=Calloc(5000,DMatrix);
    #else
    pwm=calloc(5000,sizeof(DMatrix));
    cpwm=calloc(5000,sizeof(DMatrix));
    #endif
    f =fopen(pwmfile,"rb");
    for (i=0;i<5000;i++) {
      if(readMatrix(f,&pwm[i])!=0) {break;}
      getComplementaryMotif(&pwm[i],&cpwm[i]);
    }
    Nmotif=i;
    fclose(f);

      for (j=0; j<Nmotif; j++) {
        fprintf(stdout, "\nPrint forward score contributions:\n");
        printScoreMatrix(pwm, station, trans, dx, order);
 //       fprintf(stdout, "\nPrint reverse score contributions:\n");
//        printScoreMatrix(cpwm, station, trans, dx, order);
      }
      for (j=0; j<Nmotif; j++) {
        fprintf(stdout, "\nPrint Reverse score contributions:\n");
        printScoreMatrix(cpwm, station, trans, dx, order);
 //       fprintf(stdout, "\nPrint reverse score contributions:\n");
//        printScoreMatrix(cpwm, station, trans, dx, order);
      }
      deleteBackground(station, trans);
      //deleteMatrix(&trans[i]);

    //free(station);
    //free(trans);
    for (i=0;i<Nmotif;i++) {
      deleteMatrix(&pwm[i]);
      deleteMatrix(&cpwm[i]);
    }
    #ifdef IN_R
    Free(pwm);
    Free(cpwm);
    #else
    free(pwm);
    free(cpwm);
    #endif
    //free(null.totalScore.y);
}
#endif
