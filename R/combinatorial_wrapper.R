#' Combinatrial model approximation of the number of motif hits
#'
#' This function approxmiates the distribution of the number of motif hits.
#' It sums over all combinations of obtaining k hits in a random sequence
#' of a given length using and efficient dynamic programming algorithm.
#' This function is an alternative to \code{\link{compoundPoissonDist}}.
#' However, it requires fixed-length sequences and currently
#' only supports the distribution of the number of hits when both
#' DNA strands are scanned for motif hits.
#'
#'
#' @param seqlen Integer-valued vector that defines the lengths of the
#' individual sequences. This information can be extracted for
#' given DNA sequence of interest using \code{\link{numMotifHits}}.
#' @param overlap Overlap-object that contains the overlapping 
#' hit probabilities. An Overlap-object is created by 
#' \code{\link{probOverlapHit}}
#' @return List that contains
#' \describe{
#' \item{dist}{Distribution of the number of motif hits}
#' }
#'
#' @seealso \code{\link{compoundPoissonDist}}
#' @examples
#'
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' alpha=0.001
#' motifcounterOption(alpha)
#'
#' # estimate a background model
#' bg=readBackground(seqs,1)
#'
#' # load a motif from motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute overlap probabilities
#' op=probOverlapHit(motif,bg,singlestranded=FALSE)
#' seqlen=rep(100,1)
#' # Computes the distribution of the number of motif hits
#' dist=combinatorialDist(seqlen, op)
#'
#' @export
combinatorialDist=function(seqlen, overlap) {
    overlapValid(overlap)

    # The remaining sequence must be of equal length!
    if (!all(seqlen==seqlen[1])) {
        stop("The individual sequences must be of equal
            length for dynamic programming!")
    }

    if (overlap$singlestranded==TRUE) {
        stop("Currently the combinatorial model only supports scanning
            of both DNA strands")
    }
    # Analysing too short sequences might result in biases
    # of the model
    if (seqlen[1]<70) {
        warning("The sequence length might be too short
            which could cause an offset between the true distribution
            and the combinatorial model of the number of hits")
    }
    # the maximal number of hits is given by the sequence length
    maxhits=seqlen[1]
    # allocate the distribution of the number of hits
    # across multiple sequences
    totalmaxhits=maxhits*length(seqlen)
    dist=numeric(totalmaxhits+1)
    ret=.C("motifcounter_combinatorialDist",
        overlap$alpha, overlap$beta, overlap$beta3p, overlap$beta5p,
        as.numeric(dist), as.integer(seqlen[1]),
        as.integer(maxhits), as.integer(length(seqlen)),
        length(overlap$beta),
        as.integer(overlap$singlestranded),PACKAGE="motifcounter")
    
    if (is.na(sum(ret[[5]]))) {
        # This might be the case if too stringent
        # and too long or too many sequences shall be
        # assessed for significant enrichment.
        # In this case, the combinatorial model is pushed
        # to its limits
        stop("The combinatorial model experienced numerical issues,
            please try to reduce the number or length of the DNA sequences,
            reduce the false positive probability using motifcounterOption
            or try the 'compound Poisson approximation'.")
    }
    return(list(dist=ret[[5]]))
}
