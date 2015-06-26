#include <string.h>
#ifdef IN_R
#include <R.h>
#endif
#include "forground.h"
#include "sequence.h"

#define LINELENGTH 10000

int getTableMotifWidth(FILE *f) {
  float a,c,g,t;
  int mwidth=0;
  while(fscanf(f,"%f\t%f\t%f\t%f\n",&a,&c,&g,&t)!=EOF) { mwidth++; }
  return mwidth;
}

int getTransfacMotifWidth(FILE *f) {
  int mwidth=0;
  int startreading=0;
  char line[LINELENGTH];
  while(fgets(line, LINELENGTH, f)!=NULL) { 
    if(strncmp(line,"P0",2)==0) {
      startreading=1;
      continue;
    }
    if(strncmp(line,"XX",2)==0) {
      startreading=0;
      continue;
    }
    if (startreading==1) {
      mwidth++;
    }
  }
  return mwidth;
}

int getJasparMotifWidth(FILE *f) {
  #ifdef IN_R
  error("getJasparMotifWidth ... to be implemented\n");
  #else
  fprintf(stderr, "getJasparMotifWidth ... to be implemented\n");
  #endif
  return -1;
}

void printMotif(DMatrix *pwm) {
  int i=0;
  for (i=0; i<pwm->nrow; i++) {
  #ifdef IN_R
    Rprintf("%f\t%f\t%f\t%f\n", 
    pwm->data[i*ALPHABETSIZE],
    pwm->data[i*ALPHABETSIZE+1],
    pwm->data[i*ALPHABETSIZE+2],
    pwm->data[i*ALPHABETSIZE+3]);
    #else
    printf("%f\t%f\t%f\t%f\n", 
    pwm->data[i*ALPHABETSIZE],
    pwm->data[i*ALPHABETSIZE+1],
    pwm->data[i*ALPHABETSIZE+2],
    pwm->data[i*ALPHABETSIZE+3]);
    #endif
  }
}

void getComplementaryMotif(DMatrix *pwm, DMatrix *cpwm) {
  int i=0, j=0;
  cpwm->nrow=pwm->nrow;
  cpwm->ncol=pwm->ncol;
  #ifdef IN_R
  cpwm->data=Calloc(pwm->nrow*pwm->ncol,double);
  #else
  cpwm->data=calloc(pwm->nrow*pwm->ncol,sizeof(double));
  #endif
  for (i=0; i<pwm->ncol; i++) {
    for (j=0; j<pwm->nrow; j++) {
      cpwm->data[(pwm->nrow-j-1)*ALPHABETSIZE +getComplementFromIndex(i)]=
      pwm->data[j*ALPHABETSIZE +i];
    }
  }
  return;
}

float posmin(float a, float b) {
  if (a>0 && b>0) return (a<b) ? a : b;
  else if (a>0) return (a);
  else if (b>0) return (b);
  else return -1;
}

// in tab format
void getTableMotif(FILE *f, DMatrix *m, double pseudocount) {
  int i=0;
  float a,c,g,t, min=1, n;
  while(fscanf(f,"%f\t%f\t%f\t%f\n",&a, &c,&g,&t)!=EOF) { 
   n=a+c+g+t+4*pseudocount;
   m->data[i*4]=((double)a + pseudocount)/(n);
   m->data[i*4+1]=((double)c + pseudocount)/(n);
   m->data[i*4+2]=((double)g + pseudocount)/(n);
   m->data[i*4+3]=((double)t + pseudocount)/(n);
     min=(min<posmin(posmin(m->data[i*4],m->data[i*4+1]),
     posmin(m->data[i*4+2],m->data[i*4+3]))) ? min :
     posmin(posmin(m->data[i*4],m->data[i*4+1]), posmin(m->data[i*4+2],m->data[i*4+3]));
   i++;
  }

  return;
}

void getTransfacMotif(FILE *f, DMatrix *m, double pseudocount) {
  int i=0, p;
  float a,c,g,t, n;
  char letter;
  int startreading=0;
  char line[LINELENGTH];
  while(fgets(line, LINELENGTH,f)!=NULL) { 
    if(strncmp(line,"P0",2)==0) {
      startreading=1;
      continue;
    }
    if(strncmp(line,"XX",2)==0) {
      startreading=0;
      continue;
    }
    if (startreading==1) {
      sscanf(line, "%d\t%f\t%f\t%f\t%f\t%c\n",&p, &a, &c,&g,&t,&letter);
      n=a+c+g+t+4*pseudocount;
   m->data[i*4]=((double)a + pseudocount)/(n);
   m->data[i*4+1]=((double)c + pseudocount)/(n);
   m->data[i*4+2]=((double)g + pseudocount)/(n);
   m->data[i*4+3]=((double)t + pseudocount)/(n);
  //   min=(min<posmin(posmin(m->data[i*4],m->data[i*4+1]),
 //    posmin(m->data[i*4+2],m->data[i*4+3]))) ? min :
 //    posmin(posmin(m->data[i*4],m->data[i*4+1]), posmin(m->data[i*4+2],m->data[i*4+3]));
   i++;
   }
  }

  return;
}

void getJasparMotif(FILE *f, DMatrix *m, double pseudocount) {
#ifdef IN_R
 error("getJasparMotif ...  to be implemented!\n");
 #else
 fprintf(stderr, "getJasparMotif ...  to be implemented!\n");
 #endif
 return;
}

#ifndef IN_R
int makemotif(char *fmotif, char *binmotif, char *pseudocount, char *format) {
  FILE *f;
  DMatrix *pwm, *cpwm;
  int mwidth;
  double ps;

    f=fopen(fmotif,"r");
    if (strcmp(format, "tab")==0) {
      mwidth=getTableMotifWidth(f);
    } else if (strcmp(format, "transfac")==0){
      mwidth=getTransfacMotifWidth(f);
    } else if (strcmp(format, "jaspar")==0){
      mwidth=getJasparMotifWidth(f);
    } else {
      fprintf(stderr, "this format is unknown\n");
      return -1;
    }

    f=freopen(fmotif,"r",f);

#ifdef IN_R
    pwm=Calloc(1,DMatrix);
    pwm->data=Calloc(4*mwidth,double);
    cpwm=Calloc(1,DMatrix);
    #else
    pwm=calloc(1,sizeof(DMatrix));
    pwm->data=calloc(4*mwidth,sizeof(double));
    cpwm=calloc(1,sizeof(DMatrix));
    #endif
    pwm->ncol=4;pwm->nrow=mwidth;

    cpwm->ncol=4;pwm->nrow=mwidth;

    ps=atof(pseudocount);

    //f=fopen(argv[1],"r");
    fprintf(stdout,"forward\n");
    if (strcmp(format, "tab")==0) {
      getTableMotif(f, pwm, ps);
    } else if (strcmp(format, "transfac")==0){
      getTransfacMotif(f, pwm, ps);
    } else if (strcmp(format, "jaspar")==0){
      getJasparMotif(f, pwm, ps);
    } else {
    #ifdef IN_R
      Free(pwm->data);Free(pwm);
      Free(cpwm->data);Free(cpwm);
      #else
      free(pwm->data);free(pwm);
      free(cpwm->data);free(cpwm);
      #endif
      fprintf(stderr, "this format is unknown");
      return -1;
    }
    printMotif(pwm);
    getComplementaryMotif(pwm, cpwm);
    printf("reverse\n");
    printMotif(cpwm);
    fclose(f);

    f =fopen(binmotif,"ab");
    writeMatrix(f, pwm);
    #ifdef IN_R
    Free(pwm->data);
    Free(cpwm->data);
    Free(pwm);
    Free(cpwm);
    #else
    free(pwm->data);
    free(cpwm->data);
    free(pwm);
    free(cpwm);
    #endif
    fclose(f);
  return 0;
}
#endif
