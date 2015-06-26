#include <R.h>
#include "minmaxscore.h"
//#include "

extern DMatrix *Rpwm, *Rcpwm;
extern int Rorder;
extern double *Rstation, *Rtrans;
extern double Rgran;

void Rscorerange(int *scorerange) {
  ExtremalScore fescore, rescore;
  int fmins,fmaxs, rmins,rmaxs;
  int mins, maxs;
  double dx;

  if (scorerange==NULL) {
    error("scorerange is null");
    return;
  }
  if (Rpwm==NULL || Rcpwm==NULL || Rstation==NULL || Rtrans==NULL) {
    error("load forground and background before!");
    return;
  }

  dx=Rgran;
  Rprintf("range=%d, gran=%e\n",scorerange[0],dx);
  initExtremalScore(&fescore, dx, Rpwm->nrow, Rorder);
  initExtremalScore(&rescore, dx, Rcpwm->nrow, Rorder);

  loadMinMaxScores(Rpwm, Rstation, Rtrans, &fescore);
  loadMinMaxScores(Rcpwm, Rstation, Rtrans, &rescore);
  loadIntervalSize(&fescore, NULL);
  loadIntervalSize(&rescore, NULL);

  fmins=getTotalScoreLowerBound(&fescore);
  rmins=getTotalScoreLowerBound(&rescore);
  fmaxs=getTotalScoreUpperBound(&fescore);
  rmaxs=getTotalScoreUpperBound(&rescore);
  maxs=(fmaxs>rmaxs) ? fmaxs : rmaxs;
  mins=(fmins>rmins) ? fmins : rmins;

  scorerange[0]=maxs-mins+1;

  deleteExtremalScore(&fescore);
  deleteExtremalScore(&rescore);
  Rprintf("range=%d\n",scorerange[0]);
}

#ifdef WK
void Rmaxscorethreshold(int *scorerange, char *gran, char *threshold) {
  int i, Nmotif;
  DMatrix *pwm, *cpwm;
  double *station, *trans;
  int th;
  int order;
  ExtremalScore fescore;
  ExtremalScore rescore;
  FILE *f;
  int fmin,fmax, rmin,rmax;
  double dx;

    f =fopen(bgfile,"rb");
//      fread(&seql[i],sizeof(int),1,f);
    //test=calloc(power(50000,2),sizeof(double));
    //if (test==NULL) {return;};
    //scanf("%d",&fmin);
    readBackground(f,&station, &trans,&order);
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
    //maxscore.dx=0.01;

    //initScoreMetaInfo();
    dx=(double)atof(gran);
    th=(int)roundl(atof(threshold)/dx);
    fprintf(stderr, "threshold=%d\n",th);

    // actgacca
    // 01320113

    // backward:
    // 31102310
    //
    // coml:
    // tggtcagt
    // 02231023

    initExtremalScore(&fescore, dx, pwm->nrow, order);
    initExtremalScore(&rescore, dx, pwm->nrow, order);

    minScoresPerPositionForward(pwm,station, trans,fescore.minforward, &dx, order);
    minScoresPerPositionForward(cpwm,station, trans,rescore.minforward,&dx,  order);

    minScoresPerPositionBack(pwm,trans,fescore.minbackward, &dx, order);
    minScoresPerPositionBack(cpwm,trans,rescore.minbackward, &dx, order);

    maxScoresPerPositionForward(pwm,station, trans,fescore.maxforward, &dx, order);
    maxScoresPerPositionForward(cpwm,station, trans,rescore.maxforward,&dx,  order);
    maxScoresPerPositionBack(pwm,trans,fescore.maxbackward, &dx, order);
    maxScoresPerPositionBack(cpwm,trans,rescore.maxbackward, &dx, order);

    minMotifScoreBack(pwm,station,fescore.minbackward, &dx, &fmin, order);
    minMotifScoreBack(cpwm,station,rescore.minbackward, &dx, &rmin, order);
    maxMotifScoreBack(pwm,station,fescore.maxbackward, &dx, &fmax, order);
    maxMotifScoreBack(cpwm,station,rescore.maxbackward, &dx, &rmax, order);

    loadIntervalSize(&fescore, &th);
    fprintf(stdout, "Finterval=%d\n",fescore.Smax);
    loadIntervalSize(&rescore, &th);
    fprintf(stdout, "Rinterval=%d\n",rescore.Smax);
    fprintf(stdout, "Motif1:\n");
    fprintf(stdout,"minforward:\n");
    printExtremValues(fescore.minforward,  pwm->nrow,  order);
    fprintf(stdout,"maxforward:\n");
    printExtremValues(fescore.maxforward,  pwm->nrow,  order);
    fprintf(stdout,"minbackward:\n");
    printExtremValues(fescore.minbackward,  pwm->nrow,  order);
    fprintf(stdout,"maxbackward:\n");
    printExtremValues(fescore.maxbackward,  pwm->nrow,  order);
    fprintf(stdout,"intervalstart:\n");
    printExtremValues(fescore.intervalstart,  pwm->nrow,  order);
    fprintf(stdout,"intervalend:\n");
    printExtremValues(fescore.intervalend,  pwm->nrow,  order);
    fprintf(stdout, "\n\nMotif2:\n");
    printExtremValues(rescore.minforward,  pwm->nrow,  order);
    fprintf(stdout,"maxforward:\n");
    printExtremValues(rescore.maxforward,  pwm->nrow,  order);
    fprintf(stdout,"minbackward:\n");
    printExtremValues(rescore.minbackward,  pwm->nrow,  order);
    fprintf(stdout,"maxbackward:\n");
    printExtremValues(rescore.maxbackward,  pwm->nrow,  order);
    fprintf(stdout,"intervalstart:\n");
    printExtremValues(rescore.intervalstart,  pwm->nrow,  order);
    fprintf(stdout,"intervalend:\n");
    printExtremValues(rescore.intervalend,  pwm->nrow,  order);
    fprintf(stdout, "fmin=%d\nrmin=%d\nfmax=%d\nrmax=%d\n",
    fmin,rmin,fmax,rmax);

    
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
#endif
