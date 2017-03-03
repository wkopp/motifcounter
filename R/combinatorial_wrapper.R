#' Combinatrial model approximation of the number of motif hits
#'
#' This function approxmiates the distribution of the number of motif hits.
#' To this end, it sums over all combinations of obtaining k hits 
#' in a random sequence of a given length using an efficient 
#' dynamic programming algorithm.
#'
#' This function is an alternative to \code{\link{compoundPoissonDist}}
#' which requires fixed-length sequences and currently
#' only supports the computation of 
#' the distribution of the number of hits when both
#' DNA strands are scanned for motif hits.
#' 
#' @include comppoiss_wrapper.R
#'
#' @inheritParams compoundPoissonDist
#'
#' @return List containing
#' \describe{
#' \item{dist}{Distribution of the number of hits}
#' }
#' @seealso \code{\link{compoundPoissonDist}}
#' @seealso \code{\link{numMotifHits}}
#' @seealso \code{\link{probOverlapHit}}
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#' 
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Load background model
#' bg = readBackground(seqs, 1)
#'
#' # Compute overlap probabilities
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#'
#' # Use 2 sequences of length 100 bp each
#' seqlen = rep(100, 2) 
#'
#' # Computes the combinatorial distribution of the number of motif hits
#' dist = motifcounter:::combinatorialDist(seqlen, op)
#'
combinatorialDist = function(seqlen, overlap) {
    overlapValid(overlap)
    
    # Length must be at least as long as the motif
    if (seqlen[1] - length(overlap$beta) + 1 <= 0) {
        return (list(dist = 1))
    }
    
    # The remaining sequence must be of equal length!
    if (!all(seqlen == seqlen[1])) {
        stop(
            paste(strwrap("Fixed sequence length required!
            Trim the sequences to equal length if necessary."), 
            collapse = "\n")
        )
    }
    
    if (overlap$singlestranded == TRUE) {
        stop(
            paste(strwrap("The combinatorial model only supports scanning
            of both DNA strands. Set 'singlestranded = FALSE'
            or use the compound Poisson approximation."), collapse = "\n")
        )
    }
    
    # Analysing too short sequences might result in biases
    # of the model
    if (seqlen[1] < 30) {
        warning(paste(strwrap(
            "The sequence length too short.
            Be aware that this might cause biases in the analysis.
            Use longer sequences if possible."), collapse="\n")
        )
    }
    
    # the maximal number of hits is given by the sequence length
    maxhits = seqlen[1]
    
    # allocate the distribution of the number of hits
    # across multiple sequences
    totalmaxhits = maxhits * length(seqlen)
    dist = numeric(totalmaxhits + 1)
    
    ret = .C(
        "motifcounter_combinatorialDist",
        overlap$alpha,
        overlap$beta,
        overlap$beta3p,
        overlap$beta5p,
        as.numeric(dist),
        as.integer(seqlen[1]),
        as.integer(maxhits),
        as.integer(length(seqlen)),
        length(overlap$beta),
        as.integer(overlap$singlestranded),
        PACKAGE = "motifcounter"
    )
    
    if (is.na(sum(ret[[5]]))) {
        # This might be the case if too stringent
        # and too long or too many sequences shall be
        # assessed for significant enrichment.
        # In this case, the combinatorial model is pushed
        # to its limits
        stop(
            paste(strwrap("The combinatorial model experienced numerical issues,
            please try to reduce the number or length of the DNA sequences
            or try the 'compound Poisson approximation' instead."), 
            collapse = "\n")
        )
    }
    return(list(dist = ret[[5]]))
}
