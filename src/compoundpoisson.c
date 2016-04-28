#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef IN_R
#include <R.h>
#include <Rmath.h>
#endif
#include "score1d.h"
#include "score2d.h"
#include "compoundpoisson.h"
#include "forground.h"
#include "countdist.h"
#include "background.h"
#include "overlap.h"

//#include "countdist.h"
//#include "overlap.h"

#define DSTRANDED 2
double computePoissonParameter(int seqlen, int mlen, 
   int maxclump, double alpha,double *theta) {
  int i;
  double ec=0.0;
  for (i=0; i<maxclump; i++) {
    ec+=(theta[i*DSTRANDED] + theta[i*DSTRANDED+1])*(double)(i+1);
  }
//  ec=2.15;
  //Rprintf("slen=%d, mlen=%d, alpha=%e, ec=%e, lambda=%e\n", seqlen,mlen,alpha,ec, ((double)2*(seqlen-mlen+1)*(alpha))/(ec));
  return ((double)2*(seqlen-mlen+1)*(alpha))/(ec);
}

double *initTheta(int maxclump) {
  return Calloc(maxclump*DSTRANDED, double);
}

void deleteTheta(double *theta) {
  Free(theta);
}

void computeInitialClump(double *theta, double *gamma, int mlen) {
  int i;
  theta[0]=1-gamma[mlen];
  theta[1]=1-gamma[mlen];

  for (i=1; i<mlen; i++) {
    theta[0]*=(1-gamma[i])*(1-gamma[mlen+i]);
    theta[1]*=(1-gamma[i])*(1-gamma[mlen*2+i]);
  }
}

void computeExtentionFactorsPape(double *xi, double *gamma, int mlen) {
  int k, j;
  double xik;
  // xi
  for (k=1; k<mlen; k++) {
    xik=gamma[k]/(1-gamma[k]);
    xik*=(1-gamma[mlen])/(1-gamma[mlen+k]);
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j])/((1-gamma[k+j])*(1-gamma[mlen+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j]);
    }
    #ifdef DEBUG
    Rprintf( "xi%d=%f\n", k, xik);
    #endif
    xi[0]+=xik;
  }

  // xi',3'
  for (k=0; k<mlen; k++) {
    if (k==0) {
      xik=(gamma[mlen+k])/(1-gamma[mlen+k]);
    } else {
      xik=gamma[mlen+k]/(1-gamma[mlen+k]);
    }
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen*2+j])/((1-gamma[k+j])*(1-gamma[mlen+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen*2+j]);
    }
    #ifdef DEBUG
    Rprintf( "xi%d'=%f\n", k, xik);
    #endif
    xi[1]+=xik;
  }

  // xi,5'
  for (k=1; k<mlen; k++) {
    xik=gamma[mlen*2+k]/(1-gamma[mlen*2+k]);
    xik*=(1-gamma[mlen])/(1-gamma[k]);
    for (j=1;j<mlen-k;j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j])/((1-gamma[k+j])*(1-gamma[mlen*2+k+j]));
    }
    for (j=mlen-k; j<mlen; j++) {
      xik*=(1-gamma[j])*(1-gamma[mlen+j]);
    }
    #ifdef DEBUG
    Rprintf ("xi%d,5'=%f\n", k,xik);
    #endif
    xi[2]+=xik;
  }
  #ifdef DEBUG
  Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #endif
}

// Implementation based on separating the clumps
// and convolving individual poisson distributions
#ifdef for_later_use
static void computeCompoundPoissonDistributionKopp(double lambda, 
  int maxhit, int maxclump, double *theta, double *cp) {
  int i, k;
  double normalize=0.0;
  double *pin1, *pin2;
	
	pin1=Calloc(maxhit+1, double);
  if (pin1==NULL) {
  	error("Memory-allocation in computeCompoundPoissonDistributionKopp failed");
	}
	pin2=Calloc(maxhit+1, double);
  if (pin2==NULL) {
  	error("Memory-allocation in computeCompoundPoissonDistributionKopp failed");
	}

	// init
  for (k=0;k<=maxhit;k++) {
    pin1[k]=dpois((double)k,
    		lambda*(theta[0]+theta[1]), 0);
	}

	for (i=1; i<maxclump;i++) {
    for (k=0;k*(i+1)<=maxhit;k++) {
      pin2[k*(i+1)]=dpois((double)k,
    		lambda*(theta[(i)*DSTRANDED]+theta[i*DSTRANDED+1]), 0);
	  }
	  convolute(cp, pin1,pin2,maxhit);
	  memcpy(cp,pin1,(maxhit+1)*sizeof(double));
	}
  for (i=0; i<=maxhit; i++) normalize+=cp[i];
  for (i=0; i<=maxhit; i++) cp[i]/=normalize;
  
  Free(pin1);
  Free(pin2);
}
#endif

// Implementation of Kemp et al., which was also used 
// in Pape et al.
static void computeCompoundPoissonDistributionKemp(double lambda, 
  int maxhit, int maxclump, double *theta, double *cp) {
  int i, j;
  double p;
  double normalize=0.0;

  #ifdef DEBUG
  Rprintf( "lambda=%f\n",lambda);
  #endif
  cp[0]=exp(-lambda);

  for (i=1;i<maxhit;i++) {
    p=0.0;
    j=i-maxclump+1;
    j=(j>0) ? j : 0;
    for (; j<i;j++) {
      p+=(i-j)*(theta[(i-j-1)*DSTRANDED]+theta[(i-j-1)*DSTRANDED +1])*cp[j];
    }
    cp[i]=lambda/((double)i)*p;
  }
  for (i=0; i<maxhit; i++) normalize+=cp[i];
  for (i=0; i<maxhit; i++) cp[i]/=normalize;
  

}

void computeCompoundPoissonDistribution(double lambda, 
  int maxhit, int maxclump, double *theta, double *cp) {
  computeCompoundPoissonDistributionKemp(lambda,
  		maxhit,maxclump, theta,cp);
}

void computeInitialClumpKopp(double *theta, double *beta3p, 
  double *delta, double *deltap, int mlen) {
  theta[0]=delta[mlen-1];
  theta[1]=deltap[mlen-1]*(1-beta3p[0]);
}

#define DEBUG
#undef DEBUG
void computeExtentionFactorsKopp(double *xi, 
  double *delta, double *deltap, double *beta, 
  double *beta3p, double *beta5p, int mlen) {
  int k;

  xi[1]=beta3p[0];
  for (k=1;k<mlen;k++) {
    xi[0]+=beta[k];
    xi[1]+=beta3p[k];
    xi[2]+=beta5p[k];
  }
  xi[1]*=deltap[mlen-1]/(delta[mlen-1]);
  xi[2]*=(delta[mlen-1])/deltap[mlen-1];
  #ifdef DEBUG
  Rprintf("Extention factors: xi=%e, xi3'=%e, xi5'=%e\n",xi[0],xi[1],xi[2]);
  #endif

}
void computeTheta(int maxclump, double *theta, 
  double *extention, int mlen) {
  int i;
  double total=0.0;

  total=theta[0]+theta[1];
  //Rprintf("%e\n",total);
  for (i=1;i<maxclump; i++) {
    theta[i*DSTRANDED]=
        extention[0]*theta[(i-1)*DSTRANDED] + 
         extention[2]*theta[(i-1)*DSTRANDED+1];
    theta[i*DSTRANDED+1]=
        extention[1]*theta[(i-1)*DSTRANDED] + 
         extention[0]*theta[(i-1)*DSTRANDED+1];

    total+=(theta[i*DSTRANDED]+theta[i*DSTRANDED+1]);
    #ifdef DEBUG
  Rprintf("total=%e, theta=%e, theta'=%e\n",total, theta[i*DSTRANDED], 
  theta[i*DSTRANDED+1]);
  #endif
  }

  for (i=0;i<maxclump; i++) {
    theta[i*DSTRANDED]/=total;
    theta[i*DSTRANDED+1]/=total;
    //#define DEBUG
    #ifdef DEBUG
    Rprintf( "theta%d=%e, theta%d'=%e\n",i+1,
      theta[i*DSTRANDED], i+1,theta[i*DSTRANDED+1]);
    #endif
  }
}


