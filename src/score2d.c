#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits.h>
#include <float.h>
#include <math.h>
#ifdef IN_R
#include <R.h>
#include <R_ext/Utils.h>
#endif
#include "sequence.h"
#include "background.h"
#include "scorefunctions.h"
#include "score2d.h"

void initScore2d(Score2d *s, int l) {
    s->y = (double*)R_alloc((size_t)l * l, sizeof(double));
    memset(s->y, 0, (size_t)l * l * sizeof(double));
    s->end1 = 0;
    s->end2 = 0;
    s->start1 = l;
    s->start2 = l;
}

int initScoreDistribution2d (DMatrix *theta, double *bg1,
                             MotifScore2d *result, int order) {
    int i;

    result->mlen = theta->nrow;
    result->ScoreBuffer1 = (Score2d*)R_alloc((size_t)power(ALPHABETSIZE, order),
            sizeof(Score2d));
    result->tmpScore = (Score2d*)R_alloc((size_t)power(ALPHABETSIZE, order + 1),
            sizeof(Score2d));
    memset(result->ScoreBuffer1, 0, power(ALPHABETSIZE, order) * sizeof(Score2d));
    memset(result->tmpScore, 0, power(ALPHABETSIZE, order+1) * sizeof(Score2d));

    for (i = 0; i < power(ALPHABETSIZE, order); i++) {
        initScore2d(&result->ScoreBuffer1[i], result->meta.length);
    }

    for (i = 0; i < power(ALPHABETSIZE, order + 1); i++) {
        initScore2d(&result->tmpScore[i], result->meta.length);
    }
    return 0;
}

void resetScore2d(Score2d *score, ScoreMetaInfo *meta) {
    memset(score->y, 0, meta->length * meta->length * sizeof(double));
    score->start1 = meta->length;
    score->start2 = meta->length;
    score->end1 = 0;
    score->end2 = 0;
}

void resetMotifScore2d(MotifScore2d *ms, int order) {
    int i;

    for(i = 0; i < power(ALPHABETSIZE, order); i++) {
        resetScore2d(&ms->ScoreBuffer1[i], &ms->meta);
    }
    for(i = 0; i < power(ALPHABETSIZE, order + 1); i++) {
        resetScore2d(&ms->tmpScore[i], &ms->meta);
    }
}


// a+b=c
void addScore2d(Score2d *a, Score2d *b, ScoreMetaInfo *meta) {
    int f, r;
    if (b->start1 > b->end1) return;
    if (b->start2 > b->end2) return;

    a->start1 = (a->start1 < b->start1) ? a->start1 : b->start1;
    a->end1 = (a->end1 > b->end1) ? a->end1 : b->end1;
    a->start2 = (a->start2 < b->start2) ? a->start2 : b->start2;
    a->end2 = (a->end2 > b->end2) ? a->end2 : b->end2;

    for (f = a->start1; f <= a->end1; f++) {
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(a, b, meta, f) private(r)
#endif
        for (r = a->start2; r <= a->end2; r++) {
            a->y[r * meta->length + f] += b->y[r * meta->length + f];
        }
    }
}

double getMarginalProbability2d(MotifScore2d *a, int order) {
    int i;
    double prob = 0.0;
    for(i = 0; i < power(ALPHABETSIZE, order); i++) {
        prob += a->ScoreBuffer1[i].y[0];
    }
    return prob;
}

void ShiftMultiplyScoreIndex2d(Score2d *dest, Score2d *src,
                               int *fds, int *rds, double p,
                               int flowerprev, int fupperprev, int flowercur, int fuppercur,
                               int rlowerprev, int rupperprev, int rlowercur, int ruppercur, int srange) {
    int r, f;
    int frestinterval;
    int rrestinterval;
    int leftshift1 = 0, leftshift2 = 0;

#ifndef IN_R
#ifdef DEBUG
    fprintf(stdout, "dsf=%d, fp=[%d/%d] fc=[%d/%d] ",
            *fds, flowerprev, fupperprev, flowercur, fuppercur);
    fprintf(stdout, "dsr=%d rp=[%d/%d] rc=[%d/%d] \n",
            *rds, rlowerprev, rupperprev, rlowercur, ruppercur);
#endif
#endif
    if (src->start1 > src->end1) return;
    if (src->start2 > src->end2) return;
    if (p == 0.0) return;

    // lprev and lcur are used to reconstruct the
    // absolute index, since we are only
    // using relative indexes with src->start and src->end
    //
    // compute dest->start
    // if dest->start < 0, set start to 0. probability mass
    // is dropped out then, because
    // it cannot reach the threshold anymore.
    // if dest->end > intervalsize at this position,
    // set dest->end to lcur + intervalsize.
    // probability mass beyond that cannot deceed the threshold.
    // use memmove to copy the values between dest->start and dest->end.
    // add up values that cannot deceed the threshold.
    //
    // finally multiply all values in dest->y by p.
    //fprintf(stderr, "s.s=%d, s.e=%d, lp=%d, up=%d\n lc=%d,uc=%d\n",
    // src->start, src->end, lprev,uprev,lcur,ucur);

    // forward strand
    dest->start1 = src->start1 - flowercur + flowerprev + *fds;
    dest->end1 = src->end1 - flowercur + flowerprev + *fds;
    dest->start2 = src->start2 - rlowercur + rlowerprev + *rds;
    dest->end2 = src->end2 - rlowercur + rlowerprev + *rds;

    if (dest->end1 < 0 || dest->end2 < 0) {
        dest->start1 = 1;
        dest->start2 = 1;
        dest->end1 = 0;
        dest->end2 = 0;
        return;
    }

    if (dest->start1 < 0) {
        leftshift1 = -dest->start1;
        dest->start1 = 0;
    }
    if (dest->start2 < 0) {
        leftshift2 = -dest->start2;
        dest->start2 = 0;
    }
    // forward strand
    if (dest->end1 > fuppercur - flowercur) {
        frestinterval = dest->end1 - fuppercur + flowercur;
        dest->end1 = fuppercur - flowercur;
    } else {
        frestinterval = 0;
    }

    if (dest->end2 > ruppercur - rlowercur) {
        rrestinterval = dest->end2 - ruppercur + rlowercur;
        dest->end2 = ruppercur - rlowercur;
    } else {
        rrestinterval = 0;
    }

    if (dest->start1 > fuppercur - flowercur) {
        dest->start1 = dest->end1;
        frestinterval = src->end1 - src->start1;
    }
    if (dest->start2 > ruppercur - rlowercur) {
        dest->start2 = dest->end2;
        rrestinterval = src->end2 - src->start2;
    }

    if (dest->start1 > dest->end1 && frestinterval <= 0) {
        return;
    }
    if (dest->start2 > dest->end2 && rrestinterval <= 0) {
        return;
    }

#ifndef IN_R
#ifdef DEBUG
    fprintf(stdout, "src1=[%d,%d], dest1=[%d/%d] lshift=%d, rest=%d ",
            src->start1, src->end1, dest->start1, dest->end1,
            leftshift1, frestinterval);
    fprintf(stdout, "src2=[%d,%d], dest2=[%d/%d] lshift=%d, rest=%d\n",
            src->start2, src->end2, dest->start2, dest->end2,
            leftshift2, rrestinterval);
#endif
#endif

#ifdef _OPENMP
    #pragma omp parallel for default(none) \
    shared(dest, src,srange, p, frestinterval, leftshift1, leftshift2) \
    private(f,r) collapse(2)
#endif
    for (r = dest->start2; r <= dest->end2; r++) {
        for (f = dest->start1; f <= dest->end1; f++) {
            dest->y[r * srange + f] =
                            src->y[(r - dest->start2 + src->start2 + leftshift2) * srange +
                                                                                 f - dest->start1 + src->start1 + leftshift1] * p;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for default(none) \
    shared(dest, src,srange, p, frestinterval,leftshift1,leftshift2) \
    private(f,r)
#endif
    for (r = dest->start2; r <= dest->end2; r++) {
        for (f = 0; f < frestinterval; f++) {
            dest->y[r * srange + dest->end1] +=
                            src->y[(r - dest->start2 + src->start2 + leftshift2) * srange +
                                                                                 f + src->start1 + dest->end1 - dest->start1 + 1 + leftshift1] * p;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for default(none) \
    shared(dest,p,src,srange,rrestinterval, leftshift1,leftshift2) \
    private(f, r)
#endif
    for (f = dest->start1; f <= dest->end1; f++) {
        for (r = 0; r < rrestinterval; r++) {


            dest->y[dest->end2 * srange + f] += src->y[
                    (r + src->start2 + dest->end2 - dest->start2 + 1 + leftshift2) * srange +
                    f - dest->start1 + src->start1 + leftshift1 ] * p;

        }
    }
    for (f = 0; f < frestinterval; f++) {
        for (r = 0; r < rrestinterval; r++) {
            dest->y[dest->end2 * srange + dest->end1] +=
                            src->y[(r + src->start2 + dest->end2 - dest->start2 + 1 +
                                    leftshift2) * srange +
                                   src->start1 + dest->end1 - dest->start1 + 1 + f + leftshift1] * p;
        }
    }
}

#define MAXBFLEN 12
int computeSeqenByBruteForce(int mlen1, int mlen2,
                             int shift, ExtremalScore *escore1, ExtremalScore *escore2) {
    int steps = 0, psteps;
    int i, prefix, j;
    int max, corder;
    int imin, start;
    max = (MAXBFLEN < shift + mlen2) ? MAXBFLEN : shift + mlen2;
    psteps = 0;
    corder = escore1->order;
    if (corder == 0) corder++;

    for (j = corder - 1; j < mlen2; j++) {
        if (shift + j >= mlen1) {
            for (prefix = 0;
                    prefix < power(ALPHABETSIZE, escore1->order); prefix++) {
                psteps += (getScoreUpperBound(escore2, j, prefix) -
                           getScoreLowerBound(escore2, j, prefix) + 1) * ALPHABETSIZE;
            }
        } else {
            for (prefix = 0;
                    prefix < power(ALPHABETSIZE, escore1->order); prefix++) {
                psteps += (getScoreUpperBound(escore2, j, prefix) -
                           getScoreLowerBound(escore2, j, prefix) + 1) *
                          (getScoreUpperBound(escore1, shift + j, prefix) -
                           getScoreLowerBound(escore1, shift + j, prefix) + 1) * ALPHABETSIZE;
            }
        }
    }
    imin = corder - 1;

    start = (shift > corder) ? shift : corder;
    for (i = start; i < max; i++) {
        steps = power(ALPHABETSIZE, i + 1) * (i + 1) * 2;

        j = 0;
        for (j = i + 1; j < max; j++) {
            if (j >= mlen1) {
                // final 1 dimensional convolution
                for (prefix = 0;
                        prefix < power(ALPHABETSIZE, escore1->order); prefix++) {
                    steps += (getScoreUpperBound(escore2, j - shift, prefix) -
                              getScoreLowerBound(escore2, j - shift, prefix) + 1) *
                             ALPHABETSIZE;
                }
            } else {
                for (prefix = 0;
                        prefix < power(ALPHABETSIZE, escore1->order); prefix++) {
                    steps += (getScoreUpperBound(escore2, j - shift, prefix) -
                              getScoreLowerBound(escore2, j - shift, prefix) + 1) *
                             (getScoreUpperBound(escore1, j, prefix) -
                              getScoreLowerBound(escore1, j, prefix) + 1) * ALPHABETSIZE;
                }
            }
            if (steps > psteps) {
                break;
            }
        }
        if (steps < psteps) {
            imin = i;
            psteps = steps;
        }
    }
    return imin + 1;
}

#define DEBUG
#undef DEBUG
int computeScoreDistribution2DDP_init(DMatrix *pwm1, DMatrix *pwm2,
                                      double *trans, double *station, MotifScore2d *mscore,
                                      ExtremalScore *escore1, ExtremalScore *escore2,
                                      int shift, int nsubmotif,
                                      MotifScore1d *init1d, MotifScore2d *init2d, int order) {
    int i, s1, s2, k, K, x, l;
    Score1d *pinit;
    double pext;
    int corder = order, inter;
    int pos, jointindex;
    int score2[power(ALPHABETSIZE, order + 1)];

    if(corder == 0) (corder++);

    getScoresInitialIndex(pwm2->data, station, score2, &mscore->meta.dx, order);

    if (nsubmotif == shift + order) {
        // in this case use 1DDP
        for (i = 0; i < power(ALPHABETSIZE, order); i++) {
            if (order == 0) {
                K = ALPHABETSIZE;
                k = 0;
            } else {
                k = i;
                K = i + 1;
            }
            for (; k < K; k++) {
                if (shift + order > pwm1->nrow) {
                    // extention is required
                    pext = 1.0;
                    for (x = 0; x < power(ALPHABETSIZE, shift + order - pwm1->nrow); x++) {
                        jointindex = x * power(ALPHABETSIZE, pwm1->nrow - shift) +
                                     k / power(ALPHABETSIZE, shift + order - pwm1->nrow);

                        pinit = &init1d->ScoreBuffer1[(pwm1->nrow - 1) * power(ALPHABETSIZE,
                                                                       order) + jointindex];
                        if (pinit->start < pinit->end) continue;
                        pext = pinit->y[0];
                        for (l = 0; l < shift + order - pwm1->nrow; l++) {
                            jointindex *= ALPHABETSIZE;
                            inter = k - (k / power(ALPHABETSIZE, shift + order - pwm1->nrow - l)) *
                                    power(ALPHABETSIZE, shift + order - pwm1->nrow - l);
                            inter /= power(ALPHABETSIZE, shift + order - pwm1->nrow - l - 1);
                            jointindex += inter;
                            pext *= trans[jointindex];
                            jointindex -= (jointindex / power(ALPHABETSIZE, order)) *
                                          power(ALPHABETSIZE, order);
                        }
                        if (score2[k] < getScoreLowerBound(escore2, corder - 1, i))
                            continue;
                        mscore->ScoreBuffer1[i].y[0] += pext;
                        mscore->ScoreBuffer1[i].start1 = 0;
                        mscore->ScoreBuffer1[i].end1 = 0;
                        mscore->ScoreBuffer1[i].start2 = 0;
                        mscore->ScoreBuffer1[i].end2 = 0;
                    }
#ifndef IN_R
#ifdef DEBUG
                    fprintf(stdout, "post-init2DDP-extention mode: nsub=%d m=0 i=%d ",
                            nsubmotif, i);
                    printProb(&mscore->ScoreBuffer1[i], mscore->meta.length);
#endif
#endif
                } else {
                    // just use the right index of 1DDP
                    pos = shift + order - 1;
                    if (score2[k] < getScoreLowerBound(escore2, corder - 1, i)) continue;
                    // here we need to aggregate all values which for sure exceed the
                    // score threshold
                    mscore->ScoreBuffer1[i].start1 =
                                    init1d->ScoreBuffer1[pos * power(ALPHABETSIZE, order) + i].start;
                    mscore->ScoreBuffer1[i].end1 =
                                    init1d->ScoreBuffer1[pos * power(ALPHABETSIZE, order) + i].end;
                    mscore->ScoreBuffer1[i].start2 = 0;
                    mscore->ScoreBuffer1[i].end2 = 0;

#ifdef _OPENMP
                    #pragma omp parallel for default(none) \
                    shared(mscore, init1d, i, order, pos) private(s1)
#endif
                    for (s1 = init1d->ScoreBuffer1[pos * power(ALPHABETSIZE, order) + i].start;
                            s1 <= init1d->ScoreBuffer1[pos * power(ALPHABETSIZE, order) + i].end;
                            s1++) {
                        mscore->ScoreBuffer1[i].y[s1] =
                                        init1d->ScoreBuffer1[pos * power(ALPHABETSIZE, order) + i].y[s1];
                    }
#ifndef IN_R
#ifdef DEBUG
                    fprintf(stdout, "init2DDP-normal mode: nsub=%d m=0 i=%d ",
                            nsubmotif, i);
                    printProb(&mscore->ScoreBuffer1[i], mscore->meta.length);
#endif
#endif
                }
            }
        }
    } else {
        // in this case 2DBF is used.
        // It has to be computed prior to using it here.
        for (i = 0; i < power(ALPHABETSIZE, order); i++) {
            mscore->ScoreBuffer1[i].start1 = init2d->ScoreBuffer1[i].start1;
            mscore->ScoreBuffer1[i].end1 = init2d->ScoreBuffer1[i].end1;
            mscore->ScoreBuffer1[i].start2 = init2d->ScoreBuffer1[i].start2;
            mscore->ScoreBuffer1[i].end2 = init2d->ScoreBuffer1[i].end2;

#ifdef _OPENMP
            #pragma omp parallel for default(none) \
            shared(mscore, init2d, i) private(s1,s2) collapse(2)
#endif
            for (s1 = mscore->ScoreBuffer1[i].start1;
                    s1 <= mscore->ScoreBuffer1[i].end1; s1++) {
                for (s2 = mscore->ScoreBuffer1[i].start2;
                        s2 <= mscore->ScoreBuffer1[i].end2; s2++) {
                    mscore->ScoreBuffer1[i].y[s2 * mscore->meta.length + s1] =
                                    init2d->ScoreBuffer1[i].y[s2 * mscore->meta.length + s1];
                }
            }
#ifndef IN_R
#ifdef DEBUG
            fprintf(stdout, "init2DBF: nsub=%d m=0 i=%d ", nsubmotif, i);
            printProb(&mscore->ScoreBuffer1[i], mscore->meta.length);
#endif
#endif
        }
    }
#ifndef IN_R
#ifdef DEBUG
    fprintf(stdout, "\n");
#endif
#endif
    return 0;
}

#define DEBUG
#undef DEBUG
int computeScoreDistribution2DDP(DMatrix *pwm1, DMatrix *pwm2,
                                 double *trans, double *station, MotifScore2d *mscore,
                                 ExtremalScore *escore1, ExtremalScore *escore2,
                                 int shift, int nsubmotif, int order) {
    int i, m, j, ji;
    int pub1, plb1, cub1, clb1;
    int pub2, plb2, cub2, clb2;
    double pf;
    int corder = order;
    int score1[power(ALPHABETSIZE, order + 1)],
        score2[power(ALPHABETSIZE, order + 1)];

    if(corder == 0) (corder++);

    for (m = nsubmotif; m < pwm1->nrow + shift; m++) {
        for (i = 0; i < power(ALPHABETSIZE, order); i++) {
            if (m < pwm1->nrow) {
                getScoresIndex(&pwm1->data[(m)*ALPHABETSIZE],
                               &trans[i * ALPHABETSIZE],
                               score1, &mscore->meta.dx);
            } else {
                memset(score1, 0, power(ALPHABETSIZE, order + 1)*sizeof(int));
            }
            if (m == shift) {
                getScoresIndex(&pwm2->data[(m - shift)*ALPHABETSIZE],
                               station,
                               score2, &mscore->meta.dx);
            } else {
                getScoresIndex(&pwm2->data[(m - shift)*ALPHABETSIZE],
                               &trans[i * ALPHABETSIZE],
                               score2, &mscore->meta.dx);
            }

            for (j = 0; j < ALPHABETSIZE; j++) {

                resetScore2d(&mscore->tmpScore[i * ALPHABETSIZE + j],
                             &mscore->meta);
                if (order > 0) {
                    ji = i * ALPHABETSIZE + j;
                    ji -= (ji / power(ALPHABETSIZE, order)) *
                          power(ALPHABETSIZE, order);
                } else {
                    ji = 0;
                }

                if (shift == m) {
                    plb2 = 0;
                    pub2 = 0;
                    clb2 = getScoreLowerBound(escore2, 0, ji);
                    cub2 = getScoreUpperBound(escore2, 0, ji);

                } else {
                    plb2 = getScoreLowerBound(escore2, m - shift - 1, i);
                    pub2 = getScoreUpperBound(escore2, m - shift - 1, i);
                    clb2 = getScoreLowerBound(escore2, m - shift, ji);
                    cub2 = getScoreUpperBound(escore2, m - shift, ji);

                }

                if(m < pwm1->nrow) {
                    plb1 = getScoreLowerBound(escore1, m - 1, i);
                    pub1 = getScoreUpperBound(escore1, m - 1, i);
                    clb1 = getScoreLowerBound(escore1, m, ji);
                    cub1 = getScoreUpperBound(escore1, m, ji);
                    pf = pwm1->data[m * ALPHABETSIZE + j];
                } else {
                    plb1 = getScoreLowerBound(escore1, pwm1->nrow - 1, 0);
                    pub1 = getScoreUpperBound(escore1, pwm1->nrow - 1, 0);
                    clb1 = getScoreLowerBound(escore1, pwm1->nrow - 1, ji);
                    cub1 = getScoreUpperBound(escore1, pwm1->nrow - 1, ji);
                    pf = 1;
                }

#ifndef IN_R
#ifdef DEBUG
                fprintf(stdout, "ds=[%d, %d], cb1=[%d,%d], pb1=[%d,%d], "
                        "cb1=[%d,%d], pb2=[%d,%d]\n",
                        score1[j], score2[j], plb1, pub1, clb1, cub1, plb2, pub2, clb2, cub2);
                fprintf(stdout, "%d\n", score2[j]);
                fprintf(stdout, "fmp=%d, fm=%d, m=%d, shift=%d, "
                        "i=%d, j=%d, ji=%d\n", fmp, fm, m, shift, i, j, ji);
#endif
#endif

                ShiftMultiplyScoreIndex2d(&mscore->tmpScore[i * ALPHABETSIZE + j],
                                          &mscore->ScoreBuffer1[i],
                                          &score1[j],
                                          &score2[j],
                                          mscore->meta.prob(trans[i * ALPHABETSIZE + j], pf),
                                          plb1, pub1, clb1, cub1,
                                          plb2, pub2, clb2, cub2,
                                          mscore->meta.length);

#ifndef IN_R
#ifdef DEBUG
                fprintf(stderr, "here: m=%d, i=%d, j=%d, ji=%d\n", m, i, j, ji);
                fprintf(stderr, "s.i1=[%d,%d] s.i2=[%d,%d] s.y=%f\n",
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].start1,
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].end1,
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].start2,
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].end2,
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].y[0]);
                fprintf(stdout, "y1=%f yall=%f\n",
                        mscore->tmpScore[i * power(ALPHABETSIZE, order) + j].y[0],
                        mscore->ScoreBuffer1[i].y[0]);
#endif
#endif
                R_CheckUserInterrupt();
            }
        }
        for (i = 0; i < power(ALPHABETSIZE, order); i++) {
            resetScore2d(&mscore->ScoreBuffer1[i], &mscore->meta);

            for (j = 0; j < ALPHABETSIZE; j++) {
                addScore2d(&mscore->ScoreBuffer1[i],
                           &mscore->tmpScore[j * power(ALPHABETSIZE, order) + i], &mscore->meta);
            }
#ifndef IN_R
#ifdef DEBUG
            fprintf(stdout, "m=%d i=%d ", m, i);
            printProb(&mscore->ScoreBuffer1[i], mscore->meta.length);
#endif
#endif
        }
    }

    return 0;
}

#define DEBUG
#undef DEBUG
int computeScoreDistribution2DBruteForce(
                DMatrix *pwm1, DMatrix *pwm2, double *trans,
                double *station, MotifScore2d *mscore, int nsubmotif, int shift,
                ExtremalScore *escore1,
                ExtremalScore *escore2, int order) {
    int i, n, restindex, restlength, i2, l;
    int si1, si2, skip_reject = 0;
    int score1[power(ALPHABETSIZE, order + 1)],
        score2[power(ALPHABETSIZE, order + 1)];
    double p, p0 = 0.0, p1 = 0.0, ptmp;
    int prefix, cletter = 0, suffix, nextprefix;
    int corder = order, jointindex;
    if (corder == 0) corder++;

    getScoresInitialIndex(pwm1->data, station, score1, &mscore->meta.dx, order);
    getScoresInitialIndex(pwm2->data, station, score2, &mscore->meta.dx, order);

    for (i = 0; i < power(ALPHABETSIZE, nsubmotif); i++) {
        // in which array to store. for subsequent processing
        suffix = (i / power(ALPHABETSIZE, order)) * power(ALPHABETSIZE, order);
        suffix = i - suffix;

        // stationary distribution
        cletter = i / power(ALPHABETSIZE, nsubmotif - corder);

        si1 = score1[cletter];

        p = station[cletter];

        prefix = i / power(ALPHABETSIZE, nsubmotif - order);

        if (si1 < getScoreLowerBound(escore1, corder - 1, prefix)) {
            i += power(ALPHABETSIZE, nsubmotif - corder) - 1;
            p0 += p;
            continue;
        } else if (si1 > getScoreUpperBound(escore1, corder - 1, prefix)) {
            si1 = getScoreUpperBound(escore1, 0, prefix);
        }

        restindex = i - cletter * power(ALPHABETSIZE, nsubmotif - corder);

        restlength = (nsubmotif < escore1->len) ? nsubmotif : escore1->len;
        for (n = corder; n < restlength; n++) {
            cletter = restindex / power(ALPHABETSIZE, nsubmotif - n - 1);
            nextprefix = prefix * ALPHABETSIZE + cletter;
            nextprefix -= (nextprefix / power(ALPHABETSIZE, order)) *
                          power(ALPHABETSIZE, order);

            restindex -= cletter * power(ALPHABETSIZE, nsubmotif - n - 1);
            si1 += getScoreIndex(pwm1->data[n * ALPHABETSIZE + cletter],
                                 trans[prefix * ALPHABETSIZE + cletter], mscore->meta.dx);

            ptmp = trans[prefix * ALPHABETSIZE + cletter];
            p *= ptmp;

            if (si1 < getScoreLowerBound(escore1, n, nextprefix)) {
                i += power(ALPHABETSIZE, nsubmotif - n - 1) - 1;
                skip_reject = 1;
                break;
            } else if (si1 > getScoreUpperBound(escore1, n, nextprefix)) {
                si1 = getScoreUpperBound(escore1, n, nextprefix);
            }
            prefix = nextprefix;
        }
        if (skip_reject == 1) {
            skip_reject = 0;
            p0 += p;
            continue;
        }

        if (shift + order > pwm1->nrow) {
            restlength = nsubmotif - shift + corder - 1;
            if ((restlength - nsubmotif) > 0) {
                restlength = nsubmotif;
            } else {
                l = 0;
            }
            i2 = i - (i / power(ALPHABETSIZE, restlength)) *
                 power(ALPHABETSIZE, restlength);
            for (l = 0; l < shift + order - pwm1->nrow; l++) {
                jointindex = i2 / power(ALPHABETSIZE, restlength - l - corder - 1);
                ptmp = trans[jointindex];
                p *= ptmp;
                i2 = i2 - (jointindex / power(ALPHABETSIZE, corder)) *
                     power(ALPHABETSIZE, restlength - l - 1);
            }
        }
        restlength = nsubmotif - shift;
        i2 = i - (i / power(ALPHABETSIZE, restlength)) * power(ALPHABETSIZE,
                restlength);
        cletter = i2 / power(ALPHABETSIZE, restlength - corder);
        if (order > 0) {
            si2 = score2[cletter];
        } else {
            si2 = score2[cletter];
        }

        prefix = i2 / power(ALPHABETSIZE, restlength - order);

        if (si2 < getScoreLowerBound(escore2, corder - 1, prefix)) {
            l = restindex - (restindex / power(ALPHABETSIZE, restlength - n)) *
                power(ALPHABETSIZE, restlength - n);
            i += power(ALPHABETSIZE, restlength - n - 1) - 1 - l;
            p0 += p;
            continue;
        } else if (si2 > getScoreUpperBound(escore2, corder - 1, prefix)) {
            si2 = getScoreUpperBound(escore2, corder - 1, prefix);
        }

        // prefix must be adjusted correctly!
        // in case we jumped over a couple of entries
        restindex = i2 - cletter * power(ALPHABETSIZE, restlength - corder);


        for (n = corder; n < restlength; n++) {
            cletter = restindex / power(ALPHABETSIZE, restlength - n - 1);
            nextprefix = prefix * ALPHABETSIZE + cletter;
            nextprefix -= (nextprefix / power(ALPHABETSIZE, order)) *
                          power(ALPHABETSIZE, order);

            restindex -= cletter * power(ALPHABETSIZE, restlength - n - 1);
            si2 += getScoreIndex(pwm2->data[(n) * ALPHABETSIZE + cletter],
                                 trans[prefix * ALPHABETSIZE + cletter],
                                 mscore->meta.dx);

            if ((n + shift) >= ((nsubmotif < escore1->len) ?
                                nsubmotif : escore1->len)) {
                ptmp = trans[prefix * ALPHABETSIZE + cletter];
                p *= ptmp;
            }
            if (si2 < getScoreLowerBound(escore2, n, nextprefix)) {

                l = restindex - (restindex / power(ALPHABETSIZE, restlength - n)) *
                    power(ALPHABETSIZE, restlength - n);
                i += power(ALPHABETSIZE, restlength - n - 1) - 1 - l;
                skip_reject = 1;
                break;
            } else if (si2 > getScoreUpperBound(escore2, n, nextprefix)) {
                si2 = getScoreUpperBound(escore2, n, nextprefix);
            }

            prefix = nextprefix;

        }
        R_CheckUserInterrupt();

        if (skip_reject == 1) {
            skip_reject = 0;
            p0 += p;
            continue;
        }

        // it is not necessary to keep
        if (nsubmotif < escore1->len) {
            mscore->ScoreBuffer1[suffix].y[
                            (si2 - getScoreLowerBound(escore2, nsubmotif - shift - 1, suffix))*
                            mscore->meta.length
                            + si1 - getScoreLowerBound(escore1, nsubmotif - 1, suffix)] += p;
        } else {
            mscore->ScoreBuffer1[suffix].y[
                            (si2 - getScoreLowerBound(escore2, nsubmotif - shift - 1, suffix))*
                            mscore->meta.length] += p;
        }
        p1 += p;
        mscore->ScoreBuffer1[suffix].start1 = 0;
        mscore->ScoreBuffer1[suffix].start2 = 0;
        if (nsubmotif < escore1->len) {
            mscore->ScoreBuffer1[suffix].end1 =
                            getScoreUpperBound(escore1, nsubmotif - 1, suffix) -
                            getScoreLowerBound(escore1, nsubmotif - 1, suffix);
        } else {
            mscore->ScoreBuffer1[suffix].end1 = 0;
        }
        mscore->ScoreBuffer1[suffix].end2 =
                        getScoreUpperBound(escore2, nsubmotif - shift - 1, suffix) -
                        getScoreLowerBound(escore2, nsubmotif - shift - 1, suffix);

#ifndef IN_R
#ifdef DEBUG
        fprintf(stdout, "i %d: seq1 ", i);
        printSeq(i, nsubmotif);
        fprintf(stdout, " seq2 ");
        printSeq(i - (i / power(ALPHABETSIZE, restlength))*power(ALPHABETSIZE,
                 restlength), restlength);

        fprintf(stdout, " s1 %d s2 %d p %f suffix %d\n", si1, si2, p, suffix);
#endif
#endif

    }
#ifndef IN_R
#ifdef DEBUG
    fprintf(stdout, "BF-result:\n");
    for (i = 0; i < power(ALPHABETSIZE, order); i++) {
        fprintf(stdout, "m=0 i=%d ", i);
        printProb(&mscore->ScoreBuffer1[i], mscore->meta.length);
    }
#endif
#endif

    return 0;
}

#define DEBUG
#undef DEBUG
void computeConditionalOverlappingProbabilities(DMatrix *pwm1,
        DMatrix *pwm2, double *station,
        double *trans, FILE *fout, double *pvalue,
        int *inth, double *dx, double *gamma, int order) {

    MotifScore2d null, init2d;
    MotifScore1d init1d1, init1d2;
    ExtremalScore escore1, escore2;
    ExtremalScore uescore1, uescore2;
    int intervalsize, threshold, shift, submotiflen;
    int size1, size2;
    double quantile;

    // 1D score
    initExtremalScore(&uescore1, *dx, pwm1->nrow, order);
    initExtremalScore(&uescore2, *dx, pwm1->nrow, order);
    loadMinMaxScores(pwm1, station, trans, &uescore1);
    loadMinMaxScores(pwm2, station, trans, &uescore2);
    loadIntervalSize(&uescore1, NULL);
    loadIntervalSize(&uescore2, NULL);
    size1 = maxScoreIntervalSize(&uescore1);
    size2 = maxScoreIntervalSize(&uescore2);

    size1 = getTotalScoreUpperBound(&uescore1) -
            getTotalScoreLowerBound(&uescore1) + 1;
    size2 = getTotalScoreUpperBound(&uescore2) -
            getTotalScoreLowerBound(&uescore2) + 1;
    intervalsize = (size1 > size2) ? size1 : size2;

    initScoreMetaInfo(getTotalScoreLowerBound(&uescore1),
                      getTotalScoreUpperBound(&uescore1),
                      intervalsize, *dx, &init1d1.meta);

    initScoreMetaInfo(getTotalScoreLowerBound(&uescore2),
                      getTotalScoreUpperBound(&uescore2),
                      intervalsize, *dx, &init1d2.meta);

    initScoreDistribution1d(pwm1, trans, &init1d1, order);
    initScoreDistribution1d(pwm2, trans, &init1d2, order);

    computeScoreDistribution1d(pwm1, trans,
                               station, &init1d1, &uescore1, order);
    computeScoreDistribution1d(pwm2, trans,
                               station, &init1d2, &uescore2, order);

    quantile = getQuantileWithIndex1d(&init1d1,
                                      getQuantileIndex1d(&init1d1.totalScore, *pvalue));
    threshold = (int)(quantile / (*dx));

    loadIntervalSize(&uescore1, &threshold);
    loadIntervalSize(&uescore2, &threshold);
    cutScoreRangeWithThreshold(&init1d1, &uescore1, order);
    cutScoreRangeWithThreshold(&init1d2, &uescore2, order);

    initExtremalScore(&escore1, *dx, pwm1->nrow, order);
    initExtremalScore(&escore2, *dx, pwm1->nrow, order);

    // 2D score with threshold
    loadMinMaxScores(pwm1, station, trans, &escore1);
    loadMinMaxScores(pwm2, station, trans, &escore2);

    loadIntervalSize(&escore1, &threshold);
    loadIntervalSize(&escore2, &threshold);

    intervalsize = (maxScoreIntervalSize(&escore1) >
                    maxScoreIntervalSize(&escore2)) ?
                   maxScoreIntervalSize(&escore1) : maxScoreIntervalSize(&escore2);

    initScoreMetaInfo(0, intervalsize, intervalsize, *dx, &null.meta);
    initScoreMetaInfo(0, intervalsize, intervalsize, *dx, &init2d.meta);

    initScoreDistribution2d(pwm1, trans, &null, order);
    initScoreDistribution2d(pwm1, trans, &init2d, order);


    for (shift = 0; shift < pwm1->nrow; shift++) {
        submotiflen = computeSeqenByBruteForce(pwm1->nrow,
                                               pwm1->nrow, shift, &escore1, &escore1);
#ifdef DEBUG
        submotiflen = order;
#endif
        if (submotiflen > shift + order) {
            resetMotifScore2d(&init2d, order);
            computeScoreDistribution2DBruteForce(
                            pwm1, pwm1, trans, station, &init2d, submotiflen, shift,
                            &escore1,
                            &escore1,
                            order);
        } else {
            submotiflen = shift + order;
        }

        resetMotifScore2d(&null, order);
        computeScoreDistribution2DDP_init(
                        pwm1, pwm1, trans, station, &null,
                        &escore1, &escore1, shift, submotiflen,
                        &init1d1, &init2d, order);

        computeScoreDistribution2DDP(
                        pwm1, pwm1, trans, station, &null,
                        &escore1, &escore1, shift, submotiflen, order);

        gamma[shift] = getMarginalProbability2d(&null, order);
    }

    for (shift = 0; shift < pwm1->nrow; shift++) {
        submotiflen = computeSeqenByBruteForce(pwm1->nrow,
                                               pwm2->nrow, shift, &escore1, &escore2);
#ifdef DEBUG
        submotiflen = order;
#endif
        if (submotiflen > shift + order) {
            resetMotifScore2d(&init2d, order);
            computeScoreDistribution2DBruteForce(
                            pwm1, pwm2, trans, station, &init2d, submotiflen, shift,
                            &escore1,
                            &escore2,
                            order);
        } else {
            submotiflen = shift + order;
        }

        resetMotifScore2d(&null, order);
        computeScoreDistribution2DDP_init(
                        pwm1, pwm2, trans, station, &null,
                        &escore1, &escore2, shift, submotiflen,
                        &init1d1, &init2d, order);

        computeScoreDistribution2DDP(
                        pwm1, pwm2, trans, station, &null,
                        &escore1, &escore2, shift, submotiflen, order);

#ifndef IN_R
        if(fout != NULL) fprintf(fout, "%1.10f\n",
                                     getMarginalProbability2d(&null, order));
        printf("%1.10f\n", getMarginalProbability2d(&null, order));
#endif
        gamma[pwm1->nrow + shift] = getMarginalProbability2d(&null, order);
    }

    for (shift = 0; shift < pwm1->nrow; shift++) {
        submotiflen = computeSeqenByBruteForce(pwm2->nrow,
                                               pwm1->nrow, shift, &escore2, &escore1);
#ifdef DEBUG
        submotiflen = order;
#endif
        if (submotiflen > shift + order) {
            resetMotifScore2d(&init2d, order);
            computeScoreDistribution2DBruteForce(
                            pwm2, pwm1, trans, station, &init2d, submotiflen, shift,
                            &escore2,
                            &escore1,
                            order);
        } else {
            submotiflen = shift + order;
        }

        resetMotifScore2d(&null, order);
        computeScoreDistribution2DDP_init(
                        pwm2, pwm1, trans, station, &null,
                        &escore2, &escore1, shift, submotiflen,
                        &init1d2, &init2d, order);

        computeScoreDistribution2DDP(
                        pwm2, pwm1, trans, station, &null,
                        &escore2, &escore1, shift, submotiflen, order);

#ifndef IN_R
        if (fout != NULL) fprintf(fout, "%1.10f\n",
                                      getMarginalProbability2d(&null, order));
        fprintf(stdout, "%1.10f\n", getMarginalProbability2d(&null, order));
#endif
        gamma[pwm1->nrow * 2 + shift] = getMarginalProbability2d(&null, order);
    }

}

