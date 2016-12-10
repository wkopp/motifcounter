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
    if (class(seqs)=="DNAStringSet") {
        lseq=sapply(seqs, function(seq) {
            if (Biostrings::countPattern("N",seq)>0 ||
                Biostrings::countPattern("n",seq)>0) {
                return (0)
            } else {
                return(length(seq))
            }
        })
    } else if (class(seqs)=="DNAString") {
        if (Biostrings::countPattern("N",seqs)>0 ||
            Biostrings::countPattern("n",seqs)>0) {
            lseq=0
        } else {
            lseq=length(seqs)
        }
    } else {
        stop("seqs must be a DNAStringSet or a DNAString")
    }
    return(as.vector(lseq))
}
