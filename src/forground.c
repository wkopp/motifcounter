#include <string.h>
#ifdef IN_R
#include <R.h>
#endif
#include "forground.h"
#include "sequence.h"

#define LINELENGTH 10000

int getTableMotifWidth(FILE *f) {
    float a,c,g,t;
    int i, mwidth=0;
    while(1) { 
        i=fscanf(f,"%f\t%f\t%f\t%f",&a,&c,&g,&t);
        if (feof(f)) { break;}
        if (i<=0) { 
            fclose(f);
            error("The motif does not seem to be in tab format");
        }

        if (ferror(f)) {
            fclose(f);
            error("The motif does not seem to be in tab format");
        }
        mwidth++; 
    }
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
    if(ferror(f)) {
        error("The motif does not appear to be in Transfac format");
    }
    if (mwidth==0) {
        error("The motif does not appear to be in Transfac format");
    }
    return mwidth;
}

int getMemeMotifWidth(FILE *f) {
    int mwidth=0, i;
    char line[LINELENGTH];
    while(fgets(line, LINELENGTH, f)!=NULL) { 
        if(strncmp(line,"letter-",7)==0) {
            for (i=0; i<LINELENGTH; i++) {
                if (strncmp(&line[i], "w= ", 3)==0) {
                    mwidth=atoi(&line[i+3]);
                    return mwidth;
                }
            }
        }
    }
    if(ferror(f)) {
        error("The motif does not appear to be in meme format");
    }
    if (mwidth==0) {
        error("The motif does not appear to be in meme format");
    }
    return mwidth;
}

int getJasparMotifWidth(FILE *f) {
    error("getJasparMotifWidth ... to be implemented\n");
    return -1;
}

void printMotif(DMatrix *pwm) {
    int i=0;
    for (i=0; i<pwm->nrow; i++) {
        Rprintf("%f\t%f\t%f\t%f\n", 
        pwm->data[i*ALPHABETSIZE],
        pwm->data[i*ALPHABETSIZE+1],
        pwm->data[i*ALPHABETSIZE+2],
        pwm->data[i*ALPHABETSIZE+3]);
    }
}

void getComplementaryMotif(DMatrix *pwm, DMatrix *cpwm) {
    int i=0, j=0;
    cpwm->nrow=pwm->nrow;
    cpwm->ncol=pwm->ncol;
    cpwm->data=Calloc(pwm->nrow*pwm->ncol,double);
    if(cpwm->data==NULL) {
        error("Memory-allocation in getComplementaryMotif failed");
    }
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
    while(1) {
        fscanf(f,"%f\t%f\t%f\t%f",&a, &c,&g,&t);
        if (feof(f)) { break;}
        if (ferror(f)) {
            error("The motif does not seem to be in tab format");
        }
        n=a+c+g+t+4*pseudocount;
        m->data[i*4]=((double)a + pseudocount)/(n);
        m->data[i*4+1]=((double)c + pseudocount)/(n);
        m->data[i*4+2]=((double)g + pseudocount)/(n);
        m->data[i*4+3]=((double)t + pseudocount)/(n);
            min=(min<posmin(posmin(m->data[i*4],m->data[i*4+1]),
            posmin(m->data[i*4+2],m->data[i*4+3]))) ? min :
            posmin(posmin(m->data[i*4],m->data[i*4+1]), 
                    posmin(m->data[i*4+2],m->data[i*4+3]));
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
            i++;
        }
    }

    return;
}

void getMemeMotif(FILE *f, DMatrix *m, double pseudocount) {
    int i=0, p;
    float a,c,g,t, n;
    char letter;
    int startreading=0;
    char line[LINELENGTH];
    while(fgets(line, LINELENGTH,f)!=NULL) { 
        if(strncmp(line,"letter-",7)==0) {
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
            i++;
        }
    }

    return;
}

void getJasparMotif(FILE *f, DMatrix *m, double pseudocount) {
#ifdef TODO_WK
    int i=0, p;
    float a,c,g,t, n;
    char letter;
    int startreading=0;
    char line[LINELENGTH];
    char subline[LINELENGTH];
    while(fgets(line, LINELENGTH,f)!=NULL) { 
        if(strncmp(line,">MA",3)==0) {
            startreading=1;
            n=0;
            continue;
        }
        if (startreading==1) {
            pos=strtok(line,"[");
            i=0;
            if (isspace(pos[1])==0) {
                m->data[i*4+n]=(atol(&pos[2]) + pseudocount);
                pos=&pos[2]
            } else {
                m->data[i*4+n]=(atol(&pos[1]) + pseudocount);
                pos=&pos[1]
            }
            i=1;
            while(strtok(pos," ")!=NULL) {
                pos=strtok(pos," ");
                pos=&pos[1];
                m->data[i*4+n]=(atol(&pos[1]) + pseudocount);
                i++;
            }
            n++;
        }
            
        for (i=0; i<m->nrow; i++) {
            for (j=0; j<4; j++) {
                tot+=m->data[i*4+j];
            }
            for (j=0; j<4; j++) {
                m->data[i*4+1]/=(double)tot;
            }
        }
#else
    warning("getJasparMotif not yet implemented");
#endif

    return;
}

