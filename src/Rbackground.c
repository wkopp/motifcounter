#include <R.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"

double *Rstation=NULL, *Rtrans=NULL;
int Rorder;
void RdestroyBackground();

void RreloadBackground(char **file) {
  FILE *f;
  RdestroyBackground();
  f=fopen(file[0],"rb");
  if(f== NULL) {
    error("Cannot open file %s",file[0]);
  }

  readBackground (f, &Rstation, &Rtrans, &Rorder);
  fclose(f);
}

void Rstorebg(char **file) {
	FILE *f;
  f=fopen(file[0],"wb");
  if (f==NULL) {
    error("Cannot open file %s",file[0]);
  }
  if (Rstation==NULL||Rtrans==NULL) {
    error("First a background model needs to be loaded");
  }

  writeBackground (f, Rstation, Rtrans, Rorder);
}

void Rmakebg(char **infasta, int *order) {
  FILE *f;
  double *count;

  RdestroyBackground();
  f =fopen(infasta[0],"r");
  if (f==NULL) {
    error("%s not found!",infasta[0]);
    return;
  }


  if (order[0]>0) {
  #ifdef IN_R
    Rstation=Calloc(power(ALPHABETSIZE,order[0]),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    Rtrans=Calloc(power(ALPHABETSIZE, order[0]+1),double);
  #else
    Rstation=calloc(power(ALPHABETSIZE,order[0]),sizeof(double));
    count=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    Rtrans=calloc(power(ALPHABETSIZE, order[0]+1),sizeof(double));
  #endif

    getNucleotideFrequencyFromSequence(f,count, order[0]);
    #ifdef WK
    for(i=0;i<power(ALPHABETSIZE, order[0]+1); i++) {
      Rprintf("c[%d]=%e\n",i,count[i]);
    }
    #endif

    getForwardTransition(count, Rtrans, order[0]);
    #ifdef WK
    for(i=0;i<power(ALPHABETSIZE, order[0]+1); i++) {
      Rprintf("T[%d]=%e\n",i,Rtrans[i]);
    }
    #endif
    getStationaryDistribution(Rtrans, Rstation, order[0]);

  } else {
  #ifdef IN_R
    Rstation=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    Rtrans=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
  #else
    Rstation=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    Rtrans=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    count=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
  #endif
    getNucleotideFrequencyFromSequence(f,count, order[0]);
    getForwardTransition(count, Rstation, order[0]);
    getForwardTransition(count, Rtrans, order[0]);
  }
  fclose(f);
  Rorder=order[0];

#ifdef WK
  for(i=0;i<power(ALPHABETSIZE, order[0]); i++) {
    Rprintf("mu[%d]=%e\n",i,Rstation[i]);
  }
  #endif
  #ifdef IN_R
  Free(count);
  #else
  free(count);
  #endif
}

void RgetOrder(int *o) {
  o[0]=Rorder;
}
void RgetBackground(double *station, double *trans) {
  int i;
  for (i=0; i<4; i++) {
    Rprintf("%e\n", Rstation[i]);
  }
  for (i=0; i<4; i++) {
    Rprintf("%e\n", Rtrans[i]);
  }
}

void RprintBackground() {
  int i;
  int order_;

  if (Rorder>0) order_=Rorder;
  else order_=1;
  for (i=0; i<power(ALPHABETSIZE, order_); i++) {
    Rprintf("mu(i=%d)=%e\n", i, Rstation[i]);
  }
  for (i=0; i<power(ALPHABETSIZE, order_+1); i++) {
    Rprintf("T(i=%d)=%e\n", i, Rtrans[i]);
  }
}

void RdestroyBackground() {

#ifdef IN_R
  if(Rstation) Free(Rstation);
  if(Rtrans) Free(Rtrans);
#else
  if(Rstation) free(Rstation);
  if(Rtrans) free(Rtrans);
#endif
  Rstation=NULL;
  Rtrans=NULL;
  //RdeleteBeta();
  //RdeleteDelta();
}

