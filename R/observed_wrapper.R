#' Number of sequences in a given fasta file
#'
#' The function parses a fasta file and returns the number of individual
#' sequence that it contains.
#'
#'
#' @param seqs DNAStringSet
#' @return The number of sequences
#' @examples
#'
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # scan the given sequence for the motif occurances
#' nseq=numSequences(seqs)
#' nseq
#'
#' @export
numSequences=function(seqs) {
    if (class(seqs)!="DNAStringSet") {
        stop("seqs must be a DNAStringSet")
    }
    length(seqs)
}



#' Length of sequences in a given fasta file
#'
#' The function parses a fasta file and returns a vector containing the length
#' of each individual sequence contained in the fasta file. Sequences
#' containing 'N' or 'n' are skipped from the analysis and are
#' set to length zero.
#'
#'
#' @param seqs DNAStringSet
#' @return A vector containing the lengths of each individual sequence in the
#' fasta file
#' @examples
#'
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # scan the given sequence for the motif occurances
#' lseq=lenSequences(seqs)
#' lseq
#'
#' @export
lenSequences=function(seqs) {
    if (class(seqs)!="DNAStringSet") {
        stop("seqs must be a DNAStringSet")
    }
    lseq=sapply(seqs, function(seq) {
        if (Biostrings::countPattern("N",seq)>0 ||
            Biostrings::countPattern("n",seq)>0) {
            return (0)
        } else {
            return(length(seq))
        }
    })
    return(as.vector(lseq))
}
