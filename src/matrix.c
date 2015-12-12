#ifdef IN_R
#include "R.h"
#endif
#include "matrix.h"

int writeMatrix(FILE *f, DMatrix *m) {
  fwrite(&m->nrow,sizeof(int),1,f);
  fwrite(&m->ncol,sizeof(int),1,f);
  fwrite(m->data,sizeof(double),m->nrow*m->ncol,f);
  return 0;
}

void deleteMatrix(DMatrix *m) {
  if (m->data!=NULL) {
    Free(m->data);
    m->data=NULL;
  }
}

int readMatrix(FILE *f,DMatrix *m) {
  int ret=0;
  ret=fread(&m->nrow,sizeof(int),1,f);
  if (ret<1) return -1;
  ret=fread(&m->ncol,sizeof(int),1,f);
  if (ret<1) return -1;
  m->data=Calloc(m->nrow*m->ncol, double);
  if(m->data==NULL) {
  	error("Memory-allocation in readMatrix failed");
	}
  ret=fread(m->data,sizeof(double),(m->nrow)*(m->ncol),f);
  if (ret<(m->nrow)*(m->ncol)) return -1;
  return 0;
}

int power (int base, int exp) {
  if (exp>0) {
    return base *power(base, exp-1);
  } else {
    return 1;
  }
}
