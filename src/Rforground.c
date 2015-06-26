
#include <R.h>
#include <Rinternals.h>
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

void Rmotiffromfile(char **fmotif, double *pseudocount, char **format) {
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
    if (strcmp(format[0], "tab")==0) {
      mwidth=getTableMotifWidth(f);
    } else if (strcmp(format[0], "transfac")==0){
      mwidth=getTransfacMotifWidth(f);
    } else if (strcmp(format[0], "jaspar")==0){
      mwidth=getJasparMotifWidth(f);
    } else {
      error("%s: format is unknown\n", format[0]);
      return;
    }

    f=freopen(fmotif[0],"r",f);

#ifdef IN_R
    Rpwm=Calloc(1,DMatrix);
    Rpwm->data=Calloc(4*mwidth,double);
    #else
    Rpwm=calloc(1,sizeof(DMatrix));
    Rpwm->data=calloc(4*mwidth,sizeof(double));
    #endif
    Rpwm->ncol=4;Rpwm->nrow=mwidth;

#ifdef IN_R
    Rcpwm=Calloc(1,DMatrix);
    #else
    Rcpwm=calloc(1,sizeof(DMatrix));
    #endif
    //cpwm->data=calloc(4*mwidth,sizeof(double));
    Rcpwm->ncol=4;Rcpwm->nrow=mwidth;

    ps=(double)pseudocount[0];

    //f=fopen(argv[1],"r");
    if (strcmp(format[0], "tab")==0) {
      getTableMotif(f, Rpwm, ps);
      getComplementaryMotif(Rpwm, Rcpwm);
    //  Rprintf("Loaded tab-format motif\n");
    } else if (strcmp(format[0], "transfac")==0){
      getTransfacMotif(f, Rpwm, ps);
      getComplementaryMotif(Rpwm, Rcpwm);
    //  Rprintf("Loaded transfac-format motif\n");
    } else if (strcmp(format[0], "jaspar")==0){
      getJasparMotif(f, Rpwm, ps);
      getComplementaryMotif(Rpwm, Rcpwm);
    //  Rprintf("Loaded jaspar-format motif\n");
    } else {
    #ifdef IN_R
      Free(Rpwm->data);Free(Rpwm);
      Free(Rcpwm->data);Free(Rcpwm);
      #else
      free(Rpwm->data);free(Rpwm);
      free(Rcpwm->data);free(Rcpwm);
      #endif
      error("this format is unknown\n");
    }
    fclose(f);
    if (checkPositivity(Rpwm)!=0) {
      Rdestroymotif();
      error("All elements of the motif must be strictly positive!");
    }
    //printMotif(Rpwm);
    //printMotif(Rcpwm);

  return;
}

#ifdef WK
SEXP Rgetmotif() {
  int i;
  SEXP pwm;
  if (Rpwm) return R_NilValue;
  pwm = PROTECT(allocMatrix(REALSXP, LENGTH(Rpwm->nrow), LENGTH(Rpwm->ncol)));
  for (i=0; i<pwm->nrow*Rpwm->ncol; i++) {
    REAL(pwm)[i]=asReal(pwm->data[i]);
  }
  U
}
#endif

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
  #ifdef IN_R
  Rpwm=Calloc(1,DMatrix);
  Rcpwm=Calloc(1,DMatrix);
  Rpwm->data=Calloc(nrow[0]*ncol[0],double);
  #else
  Rpwm=calloc(1,sizeof(DMatrix));
  Rcpwm=calloc(1,sizeof(DMatrix));
  Rpwm->data=calloc(nrow[0]*ncol[0],sizeof(double));
  #endif
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
#ifdef IN_R
  if (Rpwm) Free(Rpwm->data);
  if (Rcpwm) Free(Rcpwm->data);
  if (Rpwm) Free(Rpwm);
  if (Rcpwm) Free(Rcpwm);
#else
  if (Rpwm) free(Rpwm->data);
  if (Rcpwm) free(Rcpwm->data);
  if (Rpwm) free(Rpwm);
  if (Rcpwm) free(Rcpwm);
#endif
  Rpwm=NULL;
  Rcpwm=NULL;
//  RdeleteBeta();
//  RdeleteDelta();
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
