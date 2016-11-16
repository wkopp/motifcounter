#include <R.h>
#include "sequence.h"

void RnumSeqs(char ** fastafile, int *numofseqs) {
    FILE *f;
    char buffer[1024*16];
    numofseqs[0]=0;
    f =fopen(fastafile[0],"r");
    if (f==NULL) {
        error("IO-Error in RnumSeqs during opening of fasta file");
    }

    while(fgets(buffer, sizeof(buffer), f)!=NULL) {
        if (buffer[0]=='>') {
            numofseqs[0]++;
        }
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
    f =fopen(fastafile[0],"r");
    if (f==NULL) {
        error("Error in RnumSeqs opening fasta file");
    }

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
                //warning("Sequence number %d contains 'n' or"
                //"  'N' and is discarded.",iseq);
                writeseq=0;
                break;
            }
            if (writeheader==1 && buffer[i]=='\n') { 
                writeheader=0; 
                writeseq=1; 
                break; 
            }
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
