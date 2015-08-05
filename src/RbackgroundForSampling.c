#include <R.h>
#include "sequence.h"
#include "overlap.h"
#include "background.h"

double *RstationForSampling=NULL, *RtransForSampling=NULL;
int RorderForSampling;
void RdestroyBackgroundForSampling();

void RreloadBackgroundForSampling(char **file) {
  FILE *f;
  RdestroyBackgroundForSampling();
  f=fopen(file[0],"rb");
  if(f== NULL) {
    error("Cannot open file %s",file[0]);
  }

  readBackground (f, &RstationForSampling, &RtransForSampling, &RorderForSampling);
  fclose(f);
}

void RstorebgForSampling(char **file) {
	FILE *f;
  f=fopen(file[0],"wb");
  if (f==NULL) {
    error("Cannot open file %s",file[0]);
  }
  if (RstationForSampling==NULL||RtransForSampling==NULL) {
    error("First a background model needs to be loaded");
  }

  writeBackground (f, RstationForSampling, RtransForSampling, RorderForSampling);
}

void RmakebgForSampling(char **infasta, int *order) {
  FILE *f;
  double *count;

  RdestroyBackgroundForSampling();
  f =fopen(infasta[0],"r");
  if (f==NULL) {
    error("%s not found!",infasta[0]);
    return;
  }


  if (order[0]>0) {
  #ifdef IN_R
    RstationForSampling=Calloc(power(ALPHABETSIZE,order[0]),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    RtransForSampling=Calloc(power(ALPHABETSIZE, order[0]+1),double);
  #else
    RstationForSampling=calloc(power(ALPHABETSIZE,order[0]),sizeof(double));
    count=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    RtransForSampling=calloc(power(ALPHABETSIZE, order[0]+1),sizeof(double));
  #endif

    getNucleotideFrequencyFromSequence(f,count, order[0]);
    #ifdef WK
    for(i=0;i<power(ALPHABETSIZE, order[0]+1); i++) {
      Rprintf("c[%d]=%e\n",i,count[i]);
    }
    #endif

    getForwardTransition(count, RtransForSampling, order[0]);
    #ifdef WK
    for(i=0;i<power(ALPHABETSIZE, order[0]+1); i++) {
      Rprintf("T[%d]=%e\n",i,RtransForSampling[i]);
    }
    #endif
    getStationaryDistribution(RtransForSampling, RstationForSampling, order[0]);

  } else {
  #ifdef IN_R
    RstationForSampling=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    RtransForSampling=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
  #else
    RstationForSampling=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    RtransForSampling=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
    count=calloc(power(ALPHABETSIZE,order[0]+1),sizeof(double));
  #endif
    getNucleotideFrequencyFromSequence(f,count, order[0]);
    getForwardTransition(count, RstationForSampling, order[0]);
    getForwardTransition(count, RtransForSampling, order[0]);
  }
  fclose(f);
  RorderForSampling=order[0];

#ifdef WK
  for(i=0;i<power(ALPHABETSIZE, order[0]); i++) {
    Rprintf("mu[%d]=%e\n",i,RstationForSampling[i]);
  }
  #endif
  #ifdef IN_R
  Free(count);
  #else
  free(count);
  #endif
}

void RgetOrderForSampling(int *o) {
  o[0]=RorderForSampling;
}
void RgetBackgroundForSampling(double *station, double *trans) {
  int i;
  for (i=0; i<4; i++) {
    Rprintf("%e\n", RstationForSampling[i]);
  }
  for (i=0; i<4; i++) {
    Rprintf("%e\n", RtransForSampling[i]);
  }
}

void RprintBackgroundForSampling() {
  int i;
  int order_;

  if (RorderForSampling>0) order_=RorderForSampling;
  else order_=1;
  for (i=0; i<power(ALPHABETSIZE, order_); i++) {
    Rprintf("mu(i=%d)=%e\n", i, RstationForSampling[i]);
  }
  for (i=0; i<power(ALPHABETSIZE, order_+1); i++) {
    Rprintf("T(i=%d)=%e\n", i, RtransForSampling[i]);
  }
}

void RdestroyBackgroundForSampling() {

#ifdef IN_R
  if(RstationForSampling) Free(RstationForSampling);
  if(RtransForSampling) Free(RtransForSampling);
#else
  if(RstationForSampling) free(RstationForSampling);
  if(RtransForSampling) free(RtransForSampling);
#endif
  RstationForSampling=NULL;
  RtransForSampling=NULL;
  //RdeleteBeta();
  //RdeleteDelta();
}

