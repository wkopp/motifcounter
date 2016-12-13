#' Length of sequences in a given fasta file
#'
#' The function returns a vector containing the lengths
#' of each sequence contained in a set of sequences. Sequences
#' containing 'N' or 'n' are skipped from the analysis and are
#' set to length zero.
#'
#' @param seqs A DNAStringSet object
#' @return A vector containing the lengths of each individual sequences
#'
#' @examples
#' 
#' # Load sequences
#' file=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(file)
#' 
#' # Retrieve sequence lengths
#' motifcounter:::lenSequences(seqs)
#'
lenSequences=function(seqs) {
    stopifnot(class(seqs)=="DNAStringSet")

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
