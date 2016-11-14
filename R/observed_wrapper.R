#' Number of sequence in a provided fasta file
#' 
#' The function parses a fasta file and returns the number of individual
#' sequence in it.
#' 
#' 
#' @param seqfile Name of the sequence file in fasta format.
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' 
#' # scan the given sequence for the motif occurances
#' nseq=numSequences(seqfile)
#' nseq
#' 
#' @export
numSequences=function(seqfile) {
  nseq=integer(1)
  res=.C("mdist_numSeqs",as.character(seqfile),as.integer(nseq),PACKAGE="mdist")
    return (res[[2]])
}



#' Length of sequences provided in fasta file
#' 
#' The function parses a fasta file and returns a vector containing the length
#' of each individual sequence contained in the fasta file. Sequences
#' containing 'N' or 'n' are skipped and set to length zero.
#' 
#' 
#' @param seqfile Name of the sequence file in fasta format.
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' 
#' # scan the given sequence for the motif occurances
#' lseq=lenSequences(seqfile)
#' lseq
#' 
#' @export
lenSequences=function(seqfile) {
  nseq=numSequences(seqfile)
  lseq=integer(nseq)
  res=.C("mdist_lenSeqs",as.character(seqfile),as.integer(nseq),as.integer(lseq),PACKAGE="mdist")
    return (res[[3]])
}



#' Number of motif hits in a given DNA sequence
#' 
#' The function scans a given DNA sequence or a set of sequences and counts the
#' number of motif hits on both, the forward and the reverse strand. It returns
#' the length of each sequence supplied in the fasta file and the number of
#' hits observed across all sequences.
#' 
#' 
#' @param seqfile Name of the sequence file in fasta format.
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' alpha=0.001
#' mdistOption(alpha)
#' 
#' readBackground(seqfile,1)
#' readMotif(pwmfile,0.01)
#' 
#' # scan the given sequence on both strands for the motif occurances
#' noc=numMotifHits(seqfile)
#' noc
#' 
#' # scan the given sequence on a single strand for the motif occurances
#' noc=numMotifHits(seqfile,singlestranded=TRUE)
#' noc
#' 
#' @export
numMotifHits=function(seqfile, singlestranded=FALSE) {
  nseq=numSequences(seqfile)
  lseq=lenSequences(seqfile)
  noh=vector(mode="integer",length=nseq)
  res=.C("mdist_numberOfHits",as.character(seqfile),as.integer(noh),
      as.integer(nseq), as.integer(lseq),as.integer(singlestranded),PACKAGE="mdist")
  return (list(nseq=nseq, lseq=lseq,
      numofhits=res[[2]]))
}

