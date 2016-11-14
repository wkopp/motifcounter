
#include <R.h>
#include <Rinternals.h>
#include<string.h>
#include "overlap.h"
#include "forground.h"

DMatrix *Rpwm=NULL, *Rcpwm=NULL;

void Rdestroymotif();

// All elements of the pwm must be positive
int checkPositivity(DMatrix *pwm) {
	int i;
	for (i=0; i<pwm->nrow*pwm->ncol; i++) {
    if (pwm->data[i]==0.0) {
      return -1;
    }
  }
  return 0;
}

void Rmotiffromfile(char **fmotif, double *pseudocount) {
  FILE *f;
  //DMatrix *pwm, *cpwm;
  int mwidth;
  double ps;

  Rdestroymotif();

  f=fopen(fmotif[0],"r");
  if (f==NULL) {
    error("%s not found!", fmotif[0]);
    return;
  }
  if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("tab")], "tab")==0) {
    mwidth=getTableMotifWidth(f);
  } else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("transfac")], "transfac")==0){
    mwidth=getTransfacMotifWidth(f);
  } else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("jaspar")], "jaspar")==0){
    mwidth=getJasparMotifWidth(f);
  //} else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("meme")], "meme")==0){
  //  mwidth=getMemeMotifWidth(f);
  } else {
    error("%s: format is unknown. \n", fmotif[0]);
    return;
  }

  f=freopen(fmotif[0],"r",f);

  Rpwm=Calloc(1,DMatrix);
  if (Rpwm==NULL) {
    error("Memory-allocation in Rmotiffromfile failed");
    return;
  }
  Rpwm->data=Calloc(4*mwidth,double);
  if (Rpwm->data==NULL) {
    error("Memory-allocation in Rmotiffromfile failed");
    return;
  }
  Rpwm->ncol=4;Rpwm->nrow=mwidth;

  Rcpwm=Calloc(1,DMatrix);
  if (Rcpwm==NULL) {
    error("Memory-allocation in Rmotiffromfile failed");
    return;
  }
  Rcpwm->ncol=4;Rcpwm->nrow=mwidth;

  ps=(double)pseudocount[0];

  if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("tab")], "tab")==0) {
    getTableMotif(f, Rpwm, ps);
    getComplementaryMotif(Rpwm, Rcpwm);
  } else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("transfac")], "transfac")==0){
    getTransfacMotif(f, Rpwm, ps);
    getComplementaryMotif(Rpwm, Rcpwm);
  } else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("jaspar")], "jaspar")==0){
    getJasparMotif(f, Rpwm, ps);
    getComplementaryMotif(Rpwm, Rcpwm);
//  } else if (strcmp(&fmotif[0][strlen(fmotif[0])-strlen("meme")], "meme")==0){
//    getMemeMotif(f, Rpwm, ps);
//    getComplementaryMotif(Rpwm, Rcpwm);
  } else {
    Free(Rpwm->data);Free(Rpwm);
    Free(Rcpwm->data);Free(Rcpwm);
    error("this format is unknown\n");
  }
  fclose(f);
  if (checkPositivity(Rpwm)!=0) {
    Rdestroymotif();
    error("All elements of the motif must be strictly positive!");
    return;
  }

  return;
}

void Rmotiflength(int *mlen) {
  if (!Rpwm) error("No motif loaded. Use read.motif");
  mlen[0]=Rpwm->nrow;
  return;
}

void Rloadmotif(double *data, int *nrow, int *ncol) {
  int i;
  if (nrow[0]!=4) {
    error("nrows must be 4, which is the number of nucleotides used!");
    return;
  }
  Rdestroymotif();
  Rpwm=Calloc(1,DMatrix);
  Rcpwm=Calloc(1,DMatrix);
  if (Rpwm==NULL||Rcpwm==NULL) {
  	error("Memory-allocation in Rloadmotif failed");
    return;
  }
  Rpwm->data=Calloc(nrow[0]*ncol[0],double);
  if (Rpwm->data==NULL) {
  	error("Memory-allocation in Rloadmotif failed");
  	return;
  }
  Rpwm->ncol=nrow[0];Rpwm->nrow=ncol[0];
  Rcpwm->ncol=nrow[0];Rcpwm->nrow=ncol[0];
  for (i=0; i<nrow[0]*ncol[0]; i++) {
    Rpwm->data[i]=data[i];
    //Rprintf("%1.2e\n",data[i]);
  }
  if (checkPositivity(Rpwm)!=0) {
  	Rdestroymotif();
  	error("All elements of the PWM must be strictly positive!");
  }
  getComplementaryMotif(Rpwm, Rcpwm);
}


void Rdestroymotif() {
  if (Rpwm && Rpwm->data) Free(Rpwm->data);
  if (Rcpwm && Rcpwm->data) Free(Rcpwm->data);
  if (Rpwm) Free(Rpwm);
  if (Rcpwm) Free(Rcpwm);
  Rpwm=NULL;
  Rcpwm=NULL;
}

SEXP fetchMotif() {
	int i;
	SEXP Rmotif;
	if (!Rpwm) return R_NilValue;

  PROTECT(Rmotif=allocMatrix(REALSXP,Rpwm->ncol,Rpwm->nrow));
  for (i=0; i<Rpwm->nrow*Rpwm->ncol; i++) {
    REAL(Rmotif)[i]=Rpwm->data[i];
  }
  UNPROTECT(1);
  return Rmotif;
}
