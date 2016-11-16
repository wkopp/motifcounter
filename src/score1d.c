#include <string.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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

double getQuantileWithIndex1d(MotifScore1d *s, int qi) {
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
    }
    return i+x;
}

double getProbWithIndex1d(MotifScore1d *s, int iquantile) {
    int i;
    double sum=0.0;

    for (i=iquantile-s->meta.xmin; i<=s->totalScore.end; i++) {
        sum+=s->totalScore.y[i];
    }
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
        for (i=0; i<=s->meta.xmax-s->meta.xmin&&i<s->meta.length; i++) {
            fprintf(f, "%e ", (double)(s->meta.xmin+i)*s->meta.dx);
        }
        fprintf(f, "\n");
    }
    for (i=0; i<=s->meta.xmax-s->meta.xmin; i++) {
        fprintf(f, "%e ", s->totalScore.y[i]);
    }
    fprintf(f, "\n");
}


void initScore1d(Score1d *s, int l) {
    s->y=Calloc(l, double);
    if (s->y==NULL) {
        error("Memory-allocation in initScore1d failed");
    }
    s->end=0;
    s->merged=0;
    s->start=l;
}

int initScoreDistribution1d (DMatrix *theta, double *bg1, 
        MotifScore1d *result, int order) {
    int i;

    initScore1d(&result->totalScore, result->meta.length+1);

    result->mlen=theta->nrow;
    result->ScoreBuffer1=Calloc(power(ALPHABETSIZE, order)*theta->nrow, Score1d);
    result->tmpScore=Calloc(power(ALPHABETSIZE, order+1), Score1d);
    if (result->ScoreBuffer1==NULL||result->tmpScore==NULL) {
        error("Memory-allocation in initScoreDistribution1d failed");
    }

    for (i=0; i < power(ALPHABETSIZE, order)*theta->nrow; i++) {
        initScore1d(&result->ScoreBuffer1[i], result->meta.length+1);
    }

    for (i=0; i < power(ALPHABETSIZE, order+1); i++) {
        initScore1d(&result->tmpScore[i], result->meta.length+1);
    }
  return 0;
}

int deleteScoreDistribution1d(MotifScore1d *m, int order) {
    int j;

    for (j=0; j < power(ALPHABETSIZE, order)*m->mlen; j++) {
        Free(m->ScoreBuffer1[j].y);
    }
    for (j=0; j < power(ALPHABETSIZE, order+1); j++) {
        Free(m->tmpScore[j].y);
    }

    Free(m->ScoreBuffer1);
    Free(m->tmpScore);

    Free(m->totalScore.y);
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
    #ifdef _OPENMP
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

static void ShiftMultiplyScoreIndex1d(Score1d *dest, Score1d *src, 
  int *ds, double p, int lprev, int uprev, int lcur, int ucur) {
    int i;
    int restinterval;

    
    if (src->start > src->end) return;
    if (p == 0.0) return;

    // lprev and lcur are used to reconstruct the absolute 
    // index, since we are only
    // using relative indexes with src->start and src->end
    //
    // compute putative dest->start
    // if dest->start < 0, set start to 0. probability 
    // mass is dropped out then, because
    // it cannot reach the threshold anymore.
    // if dest->end > intervalsize at this position, 
    // set dest->end to lcur + intervalsize.
    // probability mass beyond that cannot deceed the threshold.
    // use memmove to copy the values between dest->start and dest->end.
    // add up values that cannot deceed the threshold.
    //
    // finally multiply all values in dest->y by p.

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

    if (dest->start>dest->end) { return; }

    memmove(&dest->y[dest->start],&src->y[src->start], 
            (dest->end-dest->start+1)*sizeof(double));

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(p,dest) private(i)
    #endif
    for (i=dest->start; i<=dest->end;i++) {
        dest->y[i]*=p;
    }
    for (i=0; i<restinterval;i++) {
        dest->y[dest->end]+=src->y[src->start+dest->end-dest->start+1+i]*p;
    }
}

int computeTotalScoreDistribution1d(MotifScore1d *mscore, 
        ExtremalScore *tm, int order) {

    int i=0, k;
    int lmin=getTotalScoreLowerBound(tm);
    int *lbound=getLastScoreLowerBound(tm);
    Score1d *a,*b;
    a=&mscore->totalScore;
    for (i=0; i<power(ALPHABETSIZE,order); i++) {
        b=&mscore->ScoreBuffer1[power(ALPHABETSIZE,order)*(mscore->mlen-1) + i];

        if (b->start>b->end) return 0;
        a->start=(a->start<b->start+lbound[i]-lmin) ? 
                        a->start : b->start+lbound[i]-lmin;
        a->end=(a->end>b->end+lbound[i]-lmin) ? 
                        a->end : b->end+lbound[i]-lmin;
        #ifdef _OPENMP
        #pragma omp parallel for default(none) \
                    shared(a,b,lbound,lmin, i) private(k)
        #endif
        for (k=b->start; k<=b->end;k++) {
            a->y[k+lbound[i]-lmin]+=b->y[k];
        }
    }
    return 0;
}

#define DEBUG
#undef DEBUG
void cutScoreRangeWithThreshold(MotifScore1d *mscore, 
        ExtremalScore *tm, int order) {
    int m, i, l;
    int start, end, rest;
    if (order==0) {
        m=0;
    } else {
        m=order-1;
    }

    for (; m<tm->len;m++) {
        for (i=0; i<power(ALPHABETSIZE, order); i++) {
            if(getScoreUpperBoundUnconstrainted(tm,m,i)<
                    getScoreLowerBound(tm, m, i)) {
                memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, 
                            order) +i].y, 0,
                mscore->meta.length*sizeof(double));
                mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].start=1;
                mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].end=0;
                continue;
            }
            start=getScoreLowerBound(tm, m, i)-
                        getScoreLowerBoundUnconstrainted(tm,m,i);
            end=getScoreUpperBound(tm, m, i)-
                        getScoreLowerBoundUnconstrainted(tm,m,i);
            rest=getScoreUpperBoundUnconstrainted(tm,m,i) - 
                        getScoreUpperBound(tm, m, i);
            if (start<0) {
                start=0; end=0; rest=0;
            }

            mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].start=0;
            mscore->ScoreBuffer1[m*power(ALPHABETSIZE, order) +i].end=
            getScoreUpperBound(tm,m,i)-getScoreLowerBound(tm, m, i);

            if (start>0) {
                if (mscore->meta.length>start) {
                    memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, 
                        order) +i].y, 0, start*sizeof(double));
                } else {
                    memset(mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, 
                        order) +i].y, 0, mscore->meta.length*sizeof(double));
                }
            }
            for (l=0; l<=end-start; l++) {
                mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) +i].y[l]=
                mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, 
                        order) +i].y[l+start];
                if (start>0) {
                    mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, 
                            order) +i].y[l+start]=0;
                }
            }
            for (l=0; l<rest; l++) {
                mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) 
                    +i].y[end-start]+= mscore->ScoreBuffer1[
                         (m)*power(ALPHABETSIZE, order) +i].y[l+end+1];
                mscore->ScoreBuffer1[(m)*power(ALPHABETSIZE, order) 
                        +i].y[l+end+1]=0;
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
        mscore->totalScore.y[0]+=mscore->ScoreBuffer1[
            (tm->len-1)*power(ALPHABETSIZE, order) +i].y[0];
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
        error("Background order cannot be longer than the motif.\n");
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
            getScoreUpperBound(tm,corder-1,i)-
            getScoreLowerBound(tm, corder-1, i);

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

                resetScore1d(&mscore->tmpScore[i*ALPHABETSIZE+j], 
                        &mscore->meta);
                if (order>0) {
                    ji=i*ALPHABETSIZE + j;
                    ji-=(ji/power(ALPHABETSIZE, order))*
                        power(ALPHABETSIZE, order);
                } else {
                    ji=0;
                }

                ShiftMultiplyScoreIndex1d(&mscore->tmpScore[i*ALPHABETSIZE+j], 
                    &mscore->ScoreBuffer1[(m-1)*power(ALPHABETSIZE, order)+i], 
                    &score[j],
                    mscore->meta.prob(trans[i*ALPHABETSIZE+j], 
                        pwm->data[m*ALPHABETSIZE +j]),
                    getScoreLowerBound(tm, m-1, i),
                    getScoreUpperBound(tm, m-1, i),
                    getScoreLowerBound(tm, m, ji),
                    getScoreUpperBound(tm, m, ji));

            }
        }
        for (i=0; i<power(ALPHABETSIZE, order); i++) {
            for (j=0; j<ALPHABETSIZE; j++) {
                addScore1d(&mscore->ScoreBuffer1[m*power(ALPHABETSIZE, 
                            order) +i], 
                &mscore->tmpScore[i +j*power(ALPHABETSIZE, 
                    order)], &mscore->meta);
            }
        }
    }
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
        error("Background order cannot be longer than the motif.\n");
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

