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
  #ifdef IN_R
    Free(m->data);
    #else
    free(m->data);
    #endif
    m->data=NULL;
  }
}

int readMatrix(FILE *f,DMatrix *m) {
  int ret=0;
  ret=fread(&m->nrow,sizeof(int),1,f);
  if (ret<1) return -1;
  ret=fread(&m->ncol,sizeof(int),1,f);
  if (ret<1) return -1;
  #ifdef IN_R
  m->data=Calloc(m->nrow*m->ncol, double);
  #else
  m->data=calloc(m->nrow*m->ncol, sizeof(double));
  #endif
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
