
#ifndef IN_R
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include "sharedtypes.h"
//#include "motifscore.h"
#include "scorefunctions.h"
#include "score1d.h"
#include "score2d.h"
#include "background.h"
#include "forground.h"
#include "minmaxscore.h"
#include "simulate.h"
#include "overlap.h"
#include "posteriorcount.h"
#include "compoundpoisson.h"

#define DHSSIZE 150
//define WSIZE 550
#define NPERM 1000
#define BSIZE 100*550
#define IDLEN 200
#define RNDLEN 1000000

//#define NUMSEQ 120

//#define DBG
#define DEBUG
#undef DEBUG
//#define ALPHABETSIZE 4
#define ONEFILE

char * getArg(int argc, char *argv[], char * par) {
  int i;
  for (i=2; i<argc; i++) {
    if (strcmp(argv[i], par)==0) {
      //fprintf(stdout, "%s\t%s\n", argv[i], argv[i+1]);
      return argv[i+1];
    }
  }
  return NULL;
}

int main(int argc, char *argv[])  
{

  srand(time(NULL));
  if (strcmp(argv[1],"makebg")==0) {
    makebg(getArg(argc, argv,"-fasta"),getArg(argc, argv,"-out"), 
    NULL, getArg(argc, argv,"-order"));
  } else if (strcmp(argv[1],"makefg")==0) {
    makemotif(getArg(argc,argv,"-pwm"),getArg(argc,argv,"-out"), getArg(argc, argv, "-pseudo"), 
    getArg(argc,argv,"-format"));
  } else if (strcmp(argv[1],"minmaxtest")==0) {
    testmax(getArg(argc,argv,"-pwm"), getArg(argc,argv,"-bg"), getArg(argc,argv,"-out"),
    getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"minmaxtestthreshold")==0) {
    testmaxthreshold(getArg(argc,argv,"-pwm"), getArg(argc,argv,"-bg"), getArg(argc,argv,"-out"),
    getArg(argc,argv,"-gran"), getArg(argc,argv, "-threshold"));
  } else if (strcmp(argv[1],"printscore")==0) {
    printscore(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"),getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"scoredist")==0) {
    scoredist(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-tout"), getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-th"),getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"scoredistbf")==0) {
    scoredist_bruteforce(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-tout"), getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-threshold"),
       getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"scoredist2d")==0) {
    scoredist2d(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-tout"), getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-gran"),
       getArg(argc,argv,"-threshold"));
  } else if (strcmp(argv[1],"scoredist2d_bfinit")==0) {
    scoredist2d_bfinit(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-tout"), getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-gran"),
       getArg(argc,argv,"-threshold"));
  } else if (strcmp(argv[1],"scoredist2dbf")==0) {
    scoredist2d_bruteforce(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-tout"), getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-gran"),
       getArg(argc,argv,"-threshold"),getArg(argc,argv,"-pvalue"));
  } else if (strcmp(argv[1],"comppoispape")==0) {
    compoundpoisson_correctedpape(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-slen"),getArg(argc,argv,"-maxhit"),
       getArg(argc,argv,"-maxclump"),
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-threshold"),
       getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"comppois")==0) {
    compoundpoisson_kopp(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-slen"),getArg(argc,argv,"-maxhit"),
       getArg(argc,argv,"-maxclump"),
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-threshold"),
       getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"posterior")==0) {
    testPosteriorProbability(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"), 
       getArg(argc,argv,"-dout"), 
       getArg(argc,argv,"-slen"),getArg(argc,argv,"-maxhit"),
       getArg(argc,argv,"-pvalue"),getArg(argc,argv,"-threshold"),
       getArg(argc,argv,"-numofseq"),
       getArg(argc,argv,"-gran"));
  } else if (strcmp(argv[1],"simcount")==0) {
    simulateCountDistribution(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"),
       getArg(argc,argv,"-dout"),
       getArg(argc,argv,"-gran"), getArg(argc,argv,"-pvalue"),
       getArg(argc,argv,"-nperm"), getArg(argc,argv,"-slen"),
       getArg(argc,argv,"-maxhit"), getArg(argc,argv, "-numofseq"));
#ifdef WKOPP
  } else if (strcmp(argv[1],"cntmotifs")==0) {
    cntmotifs(argv[2], argv[3], argv[4], argv[5], argv[6]);
  } else if (strcmp(argv[1],"cntparallel")==0) {
    cntparallel(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
#endif
  } else if (strcmp(argv[1],"simscore")==0) {
   simulateScores(getArg(argc,argv,"-bg"), getArg(argc,argv,"-pwm"),
   getArg(argc,argv,"-dout"),getArg(argc,argv,"-slen"), getArg(argc,argv,"-nperm"),
   getArg(argc,argv,"-gran"));
  } else {
    fprintf (stdout,"wrong arguments");
  }

  return 0;  
}  
#endif
