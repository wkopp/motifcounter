
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#ifdef IN_R
#include <R.h>
#endif
#include <math.h>
#include "scorefunctions.h"
#include "score1d.h"
#include "sequence.h"
#include "forground.h"
#include "background.h"
#include "minmaxscore.h"

int min(int a, int b) {
    return (a<b) ? a: b;
}

int max(int a, int b) {
    return (a>b) ? a: b;
}

int getMax(int *v, int N) {
    int i;
    int m=INT_MIN;
    for (i=0; i<N;i++) {
        if(m<v[i]) {
            m=v[i];
        }
    }
    return m;
}

int getMin(int *v, int N) {
    int i;
    int m=INT_MAX;
    for (i=0; i<N;i++) {
        if(m>v[i]) {
            m=v[i];
        }
    }
    return m;
}

int getExtrem(int *v, int N, int max) {
    if (max==1) {
        return getMax(v,N);
    } else {
        return getMin(v,N);
    }
}

int whichIsMin(int *v, int N) {
    int mi=-1,i;
    int m=INT_MAX;
    for (i=0; i<N;i++) {
        if(m>v[i]) {
            m=v[i];
            mi=i;
        }
    }
    return mi;
}

int whichIsMax(int *v, int N) {
    int mi=-1,i;
    int m=INT_MIN;
    for (i=0; i<N;i++) {
        if(m<v[i]) {
            m=v[i];
            mi=i;
        }
    }
    return mi;
}

int initExtremalScore(ExtremalScore *s, double dx, int length, int order) {
    s->dx=dx;

    s->maxforward=Calloc((length)*power(ALPHABETSIZE, order), int);
    s->maxbackward=Calloc((length)*power(ALPHABETSIZE, order), int);
    s->minforward=Calloc((length)*power(ALPHABETSIZE, order), int);
    s->minbackward=Calloc((length)*power(ALPHABETSIZE, order), int);
    s->intervalstart=Calloc((length)*power(ALPHABETSIZE, order), int);
    s->intervalend=Calloc((length)*power(ALPHABETSIZE, order), int);
    if(s->maxforward==NULL||s->maxbackward==NULL||s->minforward==NULL||
                    s->minbackward==NULL||
                    s->intervalstart==NULL||s->intervalend==NULL) {
        error("Memory-allocation in initExtremalScore failed");
    }

    s->len=length;
    s->order=order;
    return 0;
}

int deleteExtremalScore(ExtremalScore *s) {
    Free(s->maxforward);
    Free(s->maxbackward);
    Free(s->minforward);
    Free(s->minbackward);
    Free(s->intervalstart);
    Free(s->intervalend);
    return 0;
}

// computes the minimum or maximum score for a pwm, 
// the stationary distribution and a precomputed escore.
void extremMotifScoreBack(int max, DMatrix *pwm, 
  double *mono, int *escore, double *dx, int *ret, int order) {

    int m;
    int s[power(ALPHABETSIZE, order+1)];

    getScoresInitialIndex(pwm->data, mono,s, dx, order );

    if (order ==0) {
        for (m=0; m<ALPHABETSIZE; m++) {
            s[m]+=escore[0];
            s[0]=getExtrem(s,ALPHABETSIZE, max);
        }
    } else {
        for (m=0; m<power(ALPHABETSIZE, order); m++) {
            s[m]+=escore[(order-1)*power(ALPHABETSIZE, order)+m];
        }
    }
    *ret=getExtrem(s,power(ALPHABETSIZE, order), max);
}

void minMotifScoreBack(DMatrix *pwm, double *mono, 
  int *escore, double *dx, int *ret, int order) {
    extremMotifScoreBack(0, pwm, mono, escore, dx, ret, order);
    return;
}

void maxMotifScoreBack(DMatrix *pwm, double *mono, 
  int *escore, double *dx, int *ret, int order) {
    extremMotifScoreBack(1, pwm, mono, escore, dx, ret, order);
    return;
}

// computes the extremum scores per position in the backward direction.
// that is, what is the maximum/minimum score that can be acheived
// for the remaining M-k nucleotids, 
// with a certain nucleotid observed at position k-1.
// the result is stored in a 4^order x M matrix
void extremScoresPerPositionBack(int max, DMatrix *theta, double *trans, 
  int *extrema, double *dx, int order) {
    int i,m1, m2;
    int s[ALPHABETSIZE];
    int k;

    if (order>theta->nrow) {
        error("Order is cannot be higher than the motif length!\n");
        return;
    }
    for (m1=0; m1<power(ALPHABETSIZE, order); m1++) {
        extrema[power(ALPHABETSIZE, order)*(theta->nrow-1)+m1]=0;
    }

    for (i=theta->nrow-1;i>=order && i>0; i--) {
        for (m1=0; m1<power(ALPHABETSIZE, order); m1++) {
            getScoresIndex(&theta->data[ALPHABETSIZE*i], 
                &trans[m1*ALPHABETSIZE],s, dx);

            m2=m1*ALPHABETSIZE;
            m2-=(m2/power(ALPHABETSIZE, order))*power(ALPHABETSIZE, order);

            for (k=0; k<ALPHABETSIZE; k++) {
                if (order==0) {
                    s[k]+=extrema[i*power(ALPHABETSIZE, order)+m2];
                } else {
                    s[k]+=extrema[power(ALPHABETSIZE, order)*i+m2 +k];
                }
            }
            extrema[power(ALPHABETSIZE, order)*(i-1)+m1]=
                getExtrem(s,ALPHABETSIZE, max);
        }
    }
}

void maxScoresPerPositionBack(DMatrix *theta, double *trans, 
 int *result, double *dx, int order) {
   extremScoresPerPositionBack(1, theta, trans, result, dx, order);
   return;
}

void minScoresPerPositionBack(DMatrix *theta, double *trans, 
 int  *result, double *dx, int order) {
   extremScoresPerPositionBack(0, theta, trans, result, dx, order);
   return;
}

// computes the extremum scores in the forward direction
// that is, what is the maximum/minimum score after 
// k nucleotids having an certain  nucleotide at the end.
void extremScoresPerPositionForward(int max, DMatrix *theta, 
    double *station,double *trans,
    int *extrema, double *dx, int order) {
    int i,m1, m2, m, k, index;
    int *s;

    // This function needs to be motified
    // the results of the function are used to compute 
    // the actual memory requirements for the current 
    // motif, background, and chosen threshold
    if (order>theta->nrow) {
        error("Background order cannot be longer than the motif.\n");
        return;
    }
    if (order>1) {
        s=Calloc(power(ALPHABETSIZE,order), int);
    } else {
        s=Calloc(ALPHABETSIZE, int);
    }
    if(s==NULL) {
        error("Memory-allocation in extremScoresPerPositionForward failed");
    }

    getScoresInitialIndex(&theta->data[0], station,s, dx, order );

    if (order==0) {
        s[0]=getExtrem(s, ALPHABETSIZE, max);
    }
    for (m=0; m<power(ALPHABETSIZE, order); m++) {
        if (order==0) {
            extrema[m]=s[m];
        } else {
            extrema[(order-1)*power(ALPHABETSIZE, order)+m]=s[m];
        }
    }
    if (order==0) {
        i=order+1;
    } else {
        i=order;
    }

    for (;i<theta->nrow; i++) {
        for (m1=0; m1<power(ALPHABETSIZE, order); m1++) {
            index=m1/ALPHABETSIZE;
            index=m1-index*ALPHABETSIZE;
            if (order>0) {
                for (k=0; k<ALPHABETSIZE;k++) {
                    s[k]=getScoreIndex(theta->data[ALPHABETSIZE*i+index],
                        trans[m1+k*power(ALPHABETSIZE, order)], *dx);
                }
            } else {
                getScoresIndex(&theta->data[ALPHABETSIZE*i],
                trans,s,dx);
            }
            m2=m1/ALPHABETSIZE;
            for (k=0; k<ALPHABETSIZE;k++) {
                if (order==0) {
                    s[k]+=extrema[(i-1)*power(ALPHABETSIZE, order)];
                } else {
                    s[k]+=extrema[(i-1)*power(ALPHABETSIZE, order)+
                        m2+power(ALPHABETSIZE,order-1)*k];
                }
            }
            extrema[i*power(ALPHABETSIZE, order)+m1]=
                        getExtrem(s, ALPHABETSIZE, max);
        }
    }
    Free(s);
}

void maxScoresPerPositionForward(DMatrix *theta, double *station, double *trans,
 int *result, double *dx, int order) {
    extremScoresPerPositionForward(1, theta, station, trans, result, dx, order);
    return;
}

void minScoresPerPositionForward(DMatrix *theta, double *station, double *trans,
 int *result, double *dx, int order) {
    extremScoresPerPositionForward(0, theta, station, trans, result, dx, order);
    return;
}


// this function computes the interval boundaries for each position
// with the nucleotid X observed. it takes the chosen threshold into account
// thus, the intervals are shorter if for a given position and nucleotide, the 
// threshold cannot be exceeded or if it cannot be deceeded.
// at the last position, the score interval boils down to a single value.
void loadIntervalSize(ExtremalScore *escore, int *threshold) {
    int i, j;
    int lower, upper;
    int order=escore->order;
    escore->Smax=INT_MIN;
    if (order>1) {
        j=order-1;
    } else {
        j=0;
    }
    for (; j<escore->len; j++) {
        for (i=0; i<power(ALPHABETSIZE, order); i++) {
            lower=escore->minforward[j*power(ALPHABETSIZE, order) +i];
            if (threshold!=NULL && lower< 
                (*threshold-escore->maxbackward[j*power(ALPHABETSIZE, order) +i]
                )) {
                lower=*threshold-escore->maxbackward[j*power(ALPHABETSIZE, 
                        order) +i];
            } else if (threshold!=NULL && 
                    escore->maxbackward[j*power(ALPHABETSIZE, order) +i]==
                    escore->minbackward[j*power(ALPHABETSIZE, order) +i]) {
                lower=*threshold;
            }
            upper=escore->maxforward[j*power(ALPHABETSIZE, order) +i];
            if (threshold!=NULL && 
                    upper>(*threshold -
                    escore->minbackward[j*power(ALPHABETSIZE, order) +i])) {
                upper=*threshold-escore->minbackward[
                            j*power(ALPHABETSIZE, order) +i];
            } else if (threshold!=NULL && 
                    escore->maxbackward[j*power(ALPHABETSIZE, order) +i]==
                    escore->minbackward[j*power(ALPHABETSIZE, order) +i]) {
                upper=*threshold;
            }
            if (upper<lower) upper=lower;
                escore->intervalstart[j*power(ALPHABETSIZE, order)+i]=lower;
                escore->intervalend[j*power(ALPHABETSIZE, order)+i]=upper;
            if (escore->Smax<(upper-lower)) {
                escore->Smax=upper-lower;
            }
        }
    }
}

int getScoreLowerBound(ExtremalScore *escore, int pos, int nucindex) {
    return escore->intervalstart[pos*power(ALPHABETSIZE, 
                    escore->order) + nucindex];
}

int getScoreUpperBound(ExtremalScore *escore, int pos, int nucindex) {
    return escore->intervalend[pos*power(ALPHABETSIZE, 
                    escore->order) + nucindex];
}

int getScoreLowerBoundUnconstrainted(ExtremalScore *escore, 
        int pos, int nucindex) {
    return escore->minforward[pos*power(ALPHABETSIZE, 
                     escore->order) + nucindex];
}

int getScoreUpperBoundUnconstrainted(ExtremalScore *escore, 
        int pos, int nucindex) {
    return escore->maxforward[pos*power(ALPHABETSIZE, 
                    escore->order) + nucindex];
}

int getTotalScoreLowerBound(ExtremalScore *escore) {
    return getMin(&escore->intervalstart[(escore->len-1)*power(ALPHABETSIZE, 
                escore->order)], power(ALPHABETSIZE, escore->order));
}
int getTotalScoreUpperBound(ExtremalScore *escore) {
    return getMax(&escore->intervalend[(escore->len-1)*power(ALPHABETSIZE, 
                escore->order)], power(ALPHABETSIZE, escore->order));
}
int * getLastScoreLowerBound(ExtremalScore *escore) {
    return &escore->intervalstart[(escore->len-1)*power(ALPHABETSIZE, 
            escore->order)];
}

int * getScoreLowerBoundArray(ExtremalScore *escore, int pos) {
    return &escore->intervalstart[pos*power(ALPHABETSIZE, escore->order)];
}

int * getScoreUpperBoundArray(ExtremalScore *escore, int pos) {
    return &escore->intervalend[pos*power(ALPHABETSIZE, escore->order)];
}

int getScoreLowerBoundPos(ExtremalScore *escore, int pos) {
    return getMin(&escore->intervalstart[pos*power(ALPHABETSIZE, 
                escore->order)], power(ALPHABETSIZE, escore->order));
}

int getScoreUpperBoundPos(ExtremalScore *escore, int pos) {
    return getMax(&escore->intervalend[pos*power(ALPHABETSIZE, escore->order)],
        power(ALPHABETSIZE, escore->order));
}

int maxScoreIntervalSize(ExtremalScore *e) {
    return e->Smax;
}


void printScores(int *ass, int N, int iscore, double, double);

void printAllScores(DMatrix *theta, double *station, double *trans, 
    double *dx, int order) {
    int ass[theta->nrow];
    int s[power(ALPHABETSIZE, order)];
    int value, preindex;
    int i, k, j;
    double lf=0, lb=0;

    getScoresInitialIndex(&theta->data[0], station,s, dx, order );
    for (i=0; i<power(ALPHABETSIZE, theta->nrow); i++) {
        getAssignmentFromIndex(i, theta->nrow, ass);

        preindex=0;
        for (k=0; k<order; k++) {
            preindex+=ass[k]*power(ALPHABETSIZE, order-k-1);
            lf=log(theta->data[ALPHABETSIZE*k+ass[k]]);
        }
        value=s[preindex];
        lb=log(station[preindex]);

        for (j=order; j< theta->nrow;j++) {
            preindex=0;

            for (k=-order; k<0; k++) {
                preindex+=ass[j+k]*power(ALPHABETSIZE, -k);
            }
            lf+=log(theta->data[ALPHABETSIZE*(j)+ass[j]]);
            lb+=log(trans[preindex+ass[j]]);
            value+=getScoreIndex(theta->data[ALPHABETSIZE*(j)+ass[j]], 
                trans[preindex+ass[j]], *dx);
        }
        printScores(ass, theta->nrow, value, lf, lb);
    }
}

void printScores(int *ass, int N, int iscore, double lf, double lb) {
#ifndef IN_R
    char n[]="acgt";
    int i;
    for (i=0; i<N; i++) {
        printf("%c",n[ass[i]]);
    }
    printf(": %d\t%f\t%f\n", iscore, lf,lb);
#endif
    return;
}

void printExtremValues(int *e, int len, int order) {
#ifndef IN_R
    int i, j;
    for (i=0; i<len; i++) {
        for (j=0; j<power(ALPHABETSIZE, order); j++) {
            printf( "%d\t",e[i*power(ALPHABETSIZE, order)+j]);
        }
        printf("\n");
    }
#endif
}

void loadMinMaxScores(DMatrix *pwm, double *station, 
        double *trans, ExtremalScore *e) {
    int min, max;
    if (e->order>pwm->nrow) {
        error("Background order cannot be longer than the motif.\n");
        return;
    }
    minScoresPerPositionForward(pwm,station, 
            trans,e->minforward, &e->dx, e->order);
    minScoresPerPositionBack(pwm,trans,
            e->minbackward, &e->dx, e->order);
    maxScoresPerPositionForward(pwm,
            station, trans,e->maxforward, &e->dx, e->order);
    maxScoresPerPositionBack(pwm,trans,e->maxbackward, &e->dx, e->order);

    minMotifScoreBack(pwm,station,e->minbackward, &e->dx, &min, e->order);
    maxMotifScoreBack(pwm,station,e->maxbackward, &e->dx, &max, e->order);

}

