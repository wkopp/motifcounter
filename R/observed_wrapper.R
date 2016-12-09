#' Number of sequences in a given fasta file
#'
#' The function returns the number of individual
#' sequence.
#'
#'
#' @param seqs DNAStringSet
#' @return The number of sequences
numSequences=function(seqs) {
    if (class(seqs)!="DNAStringSet") {
        stop("seqs must be a DNAStringSet")
    }
    length(seqs)
}



#' Length of sequences in a given fasta file
#'
#' The function returns a vector containing the lengths
#' of each individual sequence contained in seqs. Sequences
#' containing 'N' or 'n' are skipped from the analysis and are
#' set to length zero.
#'
#'
#' @param seqs DNAStringSet
#' @return A vector containing the lengths of each individual sequences
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
