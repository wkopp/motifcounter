#' Combinatrial model approximation of the number of motif hits
#' 
#' This function approxmiates the distribution of the number of motif hits.
#' It is determined by efficiently  
#' summing over all combinations of obtaining k hits in a given sequence.
#' This function can be viewed as an alternative analytical approximation
#' to \code{\link{compoundPoissonDist}}.
#' Currently, the combinatorial model approximation is only supported
#' for scanning both strands of the DNA.
#' 
#' 
#' @param seqlen Integer-valued vector that defines the lengths of the
#' individual sequences. This information can be extracted for 
#' given DNA sequence of interest using \code{\link{numMotifHits}}.
#' Note that as the combinatorial model only supports the analysis
#' of a set of sequences of equal length, each element `seqlen`
#' must be the same.
#' @param overlap List that contains the overlapping hit probabilities which
#' are determined by \code{\link{probOverlapHit}}
#' @return List that contains 
#' \describe{
#' \item{dist}{Distribution of the number of motif hits}
#' }
#'
#' @seealso \code{\link{compoundPoissonDist}}
#' @examples
#' 
#' 
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' motiffile=system.file("extdata","x31.tab", package="mdist")
#' alpha=0.001
#' mdistOption(alpha)
#' 
#' # estimate background model from seqfile
#' readBackground(seqfile,1)
#' 
#' # load motif model from motiffile
#' readMotif(motiffile, 0.01)
#' 
#' # Compute overlap probabilies
#' op=probOverlapHit(singlestranded=FALSE)
#' seqlen=rep(100,1)
#' # Computes the distribution of the number of motif hits
#' dist=combinatorialDist(seqlen, op)
#' 
#' @export
combinatorialDist=function(seqlen, overlap) {
    # Only keep sequences that are longer than zero.
    # Zero length sequences are caused by "N"s in the sequence,
    # which are discarded from the analysis
    seqlen=seqlen[seqlen>0]

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
    ret=.C("mdist_combinatorialDist", 
        overlap$alpha, overlap$beta, overlap$beta3p, overlap$beta5p,
        as.numeric(dist), as.integer(seqlen[1]),
        as.integer(maxhits), as.integer(length(seqlen)),
        as.integer(overlap$singlestranded),PACKAGE="mdist")
    return(list(dist=ret[[5]]))
}

