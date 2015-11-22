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

void Rmakebg(char **infasta, int *order, int *nseq, int *lseq) {
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

    getNucleotideFrequencyFromSequence(f,count, order[0], nseq, lseq);

    getForwardTransition(count, Rtrans, order[0]);
    getStationaryDistribution(Rtrans, Rstation, order[0]);

  } else {
    Rstation=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    Rtrans=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    count=Calloc(power(ALPHABETSIZE,order[0]+1),double);
    getNucleotideFrequencyFromSequence(f,count, order[0], nseq,lseq);
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

void RnumSeqs(char ** fastafile, int *numofseqs) {
  FILE *f;
  char buffer[1024*16];
  numofseqs[0]=0;
  f =fopen(fastafile[0],"r");

  while(fgets(buffer, sizeof(buffer), f)!=NULL) {
    //for (i=0; i<strlen(buffer); i++) {
      if (buffer[0]=='>') {
        numofseqs[0]++;
      }
    //}
  }
  if (ferror(f)) {
    error("IO-Error in RnumSeqs");
  }
  fclose(f);
} 

void RlenSeqs(char ** fastafile, int *numofseqs, int * lseq) {
  FILE *f;
  char buffer[1024*16];
  int i, iseq=0;
  int writeheader=0, writeseq=0;
  //mofseqs[0]=0;
  f =fopen(fastafile[0],"r");

  while(fgets(buffer, sizeof(buffer), f)!=NULL) {
    for (i=0; i<strlen(buffer); i++) {
      if (buffer[i]=='>') {
        lseq[iseq]=0;
        iseq++;
        writeheader=1;
        writeseq=0;
      }
      if (writeseq==1 && buffer[i]=='\n') break;
      if (writeseq==1 && isNucleotide(buffer[i])==1) lseq[iseq-1]++;
      if (writeseq==1 && isNucleotide(buffer[i])<0) {
        lseq[iseq-1]=0;
        warning("Sequence number %d contains 'n' or 'N' and is discarded.",iseq);
        writeseq=0;
        break;
      }
      if (writeheader==1 && buffer[i]=='\n') { writeheader=0; writeseq=1; break; }
    }
  }
  if (iseq!= numofseqs[0]) {
    error("RlenSeqs: Number of sequences does not match!");
  }
  if (ferror(f)) {
    error("IO-Error in RnumSeqs");
  }
  fclose(f);
} 
