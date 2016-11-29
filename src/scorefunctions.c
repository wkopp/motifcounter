#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#endif
#include "scorefunctions.h"
#include "sequence.h"

double ProbinitPWM (double b, double *f, int index, int order) {
    double g=1;
    int o, l;

    for (o=0; o<order; o++) {
        l=index/power(ALPHABETSIZE, order-o-1);
        g*=f[o*ALPHABETSIZE+l];
        index-=l*power(ALPHABETSIZE, order-o-1);
    }
    return g;
}
double ProbinitBg (double b, double *f, int ass, int order) {
    return b;
}

double ProbPWM (double b, double f) {
    return f;
}
double ProbBg (double b, double f) {
    return b;
}

void getScoresIndex(double *P,double *Q, int *score, double *dx) {
    int i;

    for (i=0;i<ALPHABETSIZE; i++) {
        score[i]=getScoreIndex(P[i],Q[i],*dx);
    }
}

void getScoresInitialIndex(double *P,double *Q, int *score,
        double *dx, int order) {
    int i, j;
    int ass[order];
    double s;

    if (order==0) {
        order++;
    }
    for (i=0;i<power(ALPHABETSIZE, order); i++) {
        s=0;
        getAssignmentFromIndex(i, order, ass);
        for (j=0; j<order; j++) {
            s+=log(P[j*ALPHABETSIZE+ass[j]]);
        }
        s-=log(Q[i]);
        score[i]=(int)roundl(s/ (*dx));
        //fprintf(stdout, "si=%d, s=%f\n",score[i], s);
    }
}

double getScore(double P,double Q) {
    return log(P/Q);
}

int getScoreIndex(double P,double Q, double dx) {

    return (int)roundl(getScore(P,Q)/dx);
}

int initScoreMetaInfo (int smin, int smax, int intervalsize,
        double dx, ScoreMetaInfo *meta) {
    meta->length=(intervalsize)+1;
    meta->dx=dx;
    meta->xmax=smax;
    meta->xmin=smin;
    if (smax-smin>meta->length) {
        error("score range length error, len=%d, xmax=%d, xmin=%d",
                    intervalsize+1,smax,smin);
    }
    meta->zero=0;

    meta->prob=&ProbBg;
    meta->probinit=&ProbinitBg;
    return 0;
}


void scoreSequence(double *station, double *trans,
  DMatrix *pwm, char *seq, int seqlen, double *scores,
  double granularity, int order) {
    int i, j;
    int s, index;
    int score[power(ALPHABETSIZE, order+1)];

    // if the sequence contains any N's, do not process the scores
    for (i=0; i<seqlen;i++) {
        if (getNucIndex(seq[i])<0) {
            return;
        }
    }

    for (i=0; i< seqlen; i++) {
        R_CheckUserInterrupt();
        index=0;

        if (order>0) {
            getScoresInitialIndex(pwm->data, station,
                    score, &granularity, order);
            index=getIndexFromAssignment(&seq[i], order);
            s=score[index];
        } else {
            s=0;
        }

        for (j=order; j<pwm->nrow; j++) {

            index=index*ALPHABETSIZE + getNucIndex(seq[i+j]);

            s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
                trans[index],granularity);

            index-=(index/power(ALPHABETSIZE,order))*power(ALPHABETSIZE,order);
        }
        scores[i]=(double)(s*granularity);
    }
}

void scoreHistogram(double *station, double *trans,
  DMatrix *pwm, char *seq, int seqlen,
 double *dist, double granularity, int smin, int order) {
    int i, j;
    int s, index;
    int score[power(ALPHABETSIZE, order+1)];

    // if the sequence contains any N's, do not process the scores
    for (i=0; i<seqlen;i++) {
        if (getNucIndex(seq[i])<0) {
            return;
        }
    }

    for (i=0; i< seqlen-pwm->nrow+1; i++) {
        R_CheckUserInterrupt();
        index=0;

        if (order>0) {
            getScoresInitialIndex(pwm->data, station,
                    score, &granularity, order);
            index=getIndexFromAssignment(&seq[i], order);
            s=score[index];
        } else {
            s=0;
        }

        for (j=order; j<pwm->nrow; j++) {

            index=index*ALPHABETSIZE + getNucIndex(seq[i+j]);

            s+=getScoreIndex(pwm->data[j*ALPHABETSIZE +getNucIndex(seq[i+j])],
                trans[index],granularity);

            index-=(index/power(ALPHABETSIZE,order))*power(ALPHABETSIZE,order);
        }
        dist[s-smin]++;
    }
}
