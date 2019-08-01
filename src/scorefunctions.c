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

double ProbinitBg (double b, double *f, int ass, int order) {
    return b;
}

double ProbBg (double b, double f) {
    return b;
}

void getScoresIndex(double *P, double *Q, int *score, double *dx) {
    int i;

    for (i = 0; i < ALPHABETSIZE; i++) {
        score[i] = getScoreIndex(P[i], Q[i], *dx);
    }
}

void getScoresInitialIndex(double *P, double *Q, int *score,
                           double *dx, int order) {
    int i, j;
    int ass[order];
    double s;

    if (order == 0) {
        order++;
    }
    for (i = 0; i < power(ALPHABETSIZE, order); i++) {
        s = 0;
        getAssignmentFromIndex(i, order, ass);
        for (j = 0; j < order; j++) {
            s += log(P[j * ALPHABETSIZE + ass[j]]);
        }
        s -= log(Q[i]);
        score[i] = (int)roundl(s / (*dx));
        //fprintf(stdout, "si=%d, s=%f\n",score[i], s);
    }
}

double getScore(double P, double Q) {
    return log(P / Q);
}

int getScoreIndex(double P, double Q, double dx) {

    return (int)roundl(getScore(P, Q) / dx);
}

int initScoreMetaInfo (int smin, int smax, int intervalsize,
                       double dx, ScoreMetaInfo *meta) {
    meta->length = (intervalsize) + 1;
    meta->dx = dx;
    meta->xmax = smax;
    meta->xmin = smin;
    meta->zero = 0;

    meta->prob = &ProbBg;
    meta->probinit = &ProbinitBg;
    return 0;
}

// getPositionWeights
//
// This function takes a PFM and a background model and determines
// the per position log-likelihood values discretized the integer
// representation.
// Therefore the result is an integer matrix of dimensions
// motif length times alphabetsize.
//
// The purpose of this function is to compute the scores once, which
// includes using the log function. Afterwards, when scanning a sequence.
// the derived scores can be looked up in the matrix.
void getPositionWeights(double *station, double *trans, DMatrix *pfm, IMatrix *pwm,
                        double granularity, int order) {
    int j, i, index, ds;
    int initscore[power(ALPHABETSIZE, order)];

    memset(initscore, 0, power(ALPHABETSIZE, order)*sizeof(int));

    if (order > 0) {

        // Initialize with stationary distribution
        // for higher order models
        getScoresInitialIndex(pfm->data, station,
                              initscore, &granularity, order);

    }

    for (j = 0; j < pfm->nrow - order; j++) {
        for (index = 0; index < power(ALPHABETSIZE, order + 1); index++) {
            if (j == 0) {
                i = index / ALPHABETSIZE;
                pwm->data[index] = initscore[i];
            }

            i = index % ALPHABETSIZE;
            ds =  getScoreIndex(pfm->data[(j + order) * ALPHABETSIZE + i],
                                trans[index], granularity);

            pwm->data[j*power(ALPHABETSIZE, order + 1) + index] += ds;

        }
    }
}


void hitSequence(IMatrix *pwm, const char *seq, int seqlen, double *hits,
                   double granularity, int order, double threshold, ExtremalScore *escore) {
  int i, j;
  int s, index;

  // if the sequence contains any N's, do not process the scores
  if (getSequenceLength(seq, seqlen) < 0) {
    return;
  }

  for (i = 0; i < seqlen - pwm->nrow + 1 - order; i++) {
    R_CheckUserInterrupt();
    if (hasN(&seq[i], pwm->nrow + order) > 0) {
      hits[i] = NAN;
      continue;
    }
    for (j = 0, index = 0; j < order; j++) {
      index = index * ALPHABETSIZE + getNucIndex(seq[i + j]);
    }
    for (j = 0, s = 0; j < pwm->nrow; j++) {
      index = index * ALPHABETSIZE + getNucIndex(seq[i + j + order]);

      s += pwm->data[j*power(ALPHABETSIZE, order + 1) + index];
      index -= (index / power(ALPHABETSIZE, order)) * power(ALPHABETSIZE, order);

      if ((double)(s + escore->maxbackward[(j+order)*power(ALPHABETSIZE, order) + index])*granularity < threshold) {
        hits[i] = 0;
        break;
      }
      if ((double)(s + escore->minbackward[(j+order)*power(ALPHABETSIZE, order) + index])*granularity >= threshold) {
        hits[i] = 1;
        break;
      }
    }
    if ((double)(s*granularity) >= (double)(threshold)) hits[i] = 1.;

  }
}


void matchCount(IMatrix *pwm, const char *seq, int seqlen, double *nhits,
                 double granularity, int order,
                 double threshold, ExtremalScore *escore,
                 int ignore_ns) {
  int i, j;
  int s, index;

  // if the sequence contains any N's, do not process the scores
  if (getSequenceLength(seq, seqlen) < 0) {
    return;
  }

  for (i = 0; i < seqlen - pwm->nrow + 1 - order; i++) {
    R_CheckUserInterrupt();
    if (hasN(&seq[i], pwm->nrow + order) > 0) {
      if (ignore_ns == 0) {
         nhits[0] = NAN;
         break;
      }
      continue;
    }
    for (j = 0, index = 0; j < order; j++) {
      index = index * ALPHABETSIZE + getNucIndex(seq[i + j]);
    }
    for (j = 0, s = 0; j < pwm->nrow; j++) {
      index = index * ALPHABETSIZE + getNucIndex(seq[i + j + order]);

      s += pwm->data[j*power(ALPHABETSIZE, order + 1) + index];
      index -= (index / power(ALPHABETSIZE, order)) * power(ALPHABETSIZE, order);

      if ((double)(s + escore->maxbackward[(j+order)*power(ALPHABETSIZE, order) + index])*granularity < threshold) {
        break;
      }
      if ((double)(s + escore->minbackward[(j+order)*power(ALPHABETSIZE, order) + index])*granularity >= threshold) {
        nhits[0] += 1;
        break;
      }
    }
    if ((double)(s*granularity) >= (double)(threshold)) nhits[0] = 1.;

  }
}


void scoreSequence(IMatrix *pwm, const char *seq, int seqlen, double *scores,
                   double granularity, int order) {
    int i, j;
    int s, index;

    // if the sequence contains any N's, do not process the scores
    if (getSequenceLength(seq, seqlen) < 0) {
       return;
    }

    for (i = 0; i < seqlen - pwm->nrow + 1 - order; i++) {
        R_CheckUserInterrupt();
        if (hasN(&seq[i], pwm->nrow + order) > 0) {
            scores[i] = NAN;
            continue;
        }
        for (j = 0, index = 0; j < order; j++) {
            index = index * ALPHABETSIZE + getNucIndex(seq[i + j]);
        }
        // pwm.nrow equals pfm.nrow - order
        for (j = 0, s = 0; j < pwm->nrow; j++) {
            index = index * ALPHABETSIZE + getNucIndex(seq[i + j + order]);

            s += pwm->data[j*power(ALPHABETSIZE, order + 1) + index];
            index -= (index / power(ALPHABETSIZE, order)) * power(ALPHABETSIZE, order);

        }
        scores[i] = (double)(s * granularity);
    }
}


void scoreHistogram(double *station, double *trans,
                    DMatrix *pwm, const char *seq, int seqlen,
                    double *dist, double granularity, int smin, int order) {
    int i, j;
    int s, index;
    int score[power(ALPHABETSIZE, order + 1)];

    for (i = 0; i < seqlen - pwm->nrow + 1; i++) {
        R_CheckUserInterrupt();
        index = 0;
        if (hasN(&seq[i], pwm->nrow) > 0) {
            continue;
        }


        if (order > 0) {
            getScoresInitialIndex(pwm->data, station,
                                  score, &granularity, order);
            index = getIndexFromAssignment(&seq[i], order);
            s = score[index];
        } else {
            s = 0;
        }

        for (j = order; j < pwm->nrow; j++) {

            index = index * ALPHABETSIZE + getNucIndex(seq[i + j]);

            s += getScoreIndex(pwm->data[j * ALPHABETSIZE + getNucIndex(seq[i + j])],
                               trans[index], granularity);

            index -= (index / power(ALPHABETSIZE, order)) * power(ALPHABETSIZE, order);
        }
        dist[s - smin]++;
    }
}
