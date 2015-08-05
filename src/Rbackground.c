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
  Rprintf("makebg: order=%d\n",order[0]);

  if (order[0]>0) {
    Rstation=Calloc(power(ALPHABETSIZE,order[0]),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    Rtrans=Calloc(power(ALPHABETSIZE, order[0]+1),double);

    getNucleotideFrequencyFromSequence(f,count, order[0]);

    getForwardTransition(count, Rtrans, order[0]);
    getStationaryDistribution(Rtrans, Rstation, order[0]);

  } else {
    Rstation=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    Rtrans=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
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
  Free(count);
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

  if (Rorder>0) {
    for (i=0; i<power(ALPHABETSIZE, Rorder); i++) {
      Rprintf("mu(i=%d)=%e\n", i, Rstation[i]);
    }
    for (i=0; i<power(ALPHABETSIZE, Rorder+1); i++) {
      Rprintf("T(i=%d)=%e\n", i, Rtrans[i]);
    }
	} else {

    for (i=0; i<power(ALPHABETSIZE, 1); i++) {
      Rprintf("mu(i=%d)=%e\n", i, Rstation[i]);
    }
    for (i=0; i<power(ALPHABETSIZE, 1); i++) {
      Rprintf("T(i=%d)=%e\n", i, Rtrans[i]);
    }
	}
}

void RdestroyBackground() {

  if(Rstation) Free(Rstation);
  if(Rtrans) Free(Rtrans);
  Rstation=NULL;
  Rtrans=NULL;
}

