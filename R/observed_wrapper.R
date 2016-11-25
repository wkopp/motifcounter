#' Number of sequences in a given fasta file
#' 
#' The function parses a fasta file and returns the number of individual
#' sequence that it contains.
#' 
#' 
#' @param seqfile Fasta-file name
#' @return The number of sequences contained in the fasta file
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' 
#' # scan the given sequence for the motif occurances
#' nseq=numSequences(seqfile)
#' nseq
#' 
#' @export
numSequences=function(seqfile) {
    nseq=integer(1)
    res=.C("motifcounter_numSeqs",as.character(seqfile),
        as.integer(nseq),PACKAGE="motifcounter")
    return (res[[2]])
}



#' Length of sequences in a given fasta file
#' 
#' The function parses a fasta file and returns a vector containing the length
#' of each individual sequence contained in the fasta file. Sequences
#' containing 'N' or 'n' are skipped from the analysis and are 
#' set to length zero.
#' 
#' 
#' @param seqfile Fasta-file name
#' @return A vector containing the lengths of each individual sequence in the
#' fasta file
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' 
#' # scan the given sequence for the motif occurances
#' lseq=lenSequences(seqfile)
#' lseq
#' 
#' @export
lenSequences=function(seqfile) {
    nseq=numSequences(seqfile)
    lseq=integer(nseq)
    res=.C("motifcounter_lenSeqs",as.character(seqfile),
        as.integer(nseq),as.integer(lseq),PACKAGE="motifcounter")
    return (res[[3]])
}



#' Number of motif hits in a given DNA sequence
#' 
#' This function scans the DNA sequences contained in the fasta file
#' and counts the number of motif hits using the score threshold
#' that is associated with the false positive probability 'alpha' 
#' (see \code{\link{motifcounterOption}}. The function can be used
#' to count motif hits on one or both strands, respectively.
#' 
#' 
#' @param seqfile Fasta-file name
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits
#' @return A list containing 
#' \describe{
#' \item{nseq}{Number of individual sequences}
#' \item{lseq}{Vector of individual sequence lengths}
#' \item{numofhits}{Vector of the number of hits in each individual sequence}
#' }
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' 
#' # Set false positive probability
#' alpha=0.001
#' motifcounterOption(alpha)
#' 
#' # Estimate order-1 background model
#' readBackground(seqfile,1)
#' # read PFM from file
#' readMotif(motiffile,0.01)
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
    res=.C("motifcounter_numberOfHits",as.character(seqfile),as.integer(noh),
        as.integer(nseq), as.integer(lseq),
        as.integer(singlestranded),PACKAGE="motifcounter")
    return (list(nseq=nseq, lseq=lseq,
        numofhits=res[[2]]))
}

