#' Compound Poisson Approximation
#'
#' This function computes the distribution of the number of motif hits
#' that emerges from a random DNA sequence of a given length.
#' The distribution can be determined in two alternative ways:
#' 1) We provide a re-implemented version of the algorithm that was
#' described in Pape et al. \emph{Compound poisson approximation
#' of the number of occurrences of a position
#' frequency matrix (PFM) on both strands.} 2008. 
#' The main purpose of this method concerns benchmarking an improved 
#' approximation.
#' In contrast to the original model, this implementation 
#' can be used with general order-d Markov models.
#' To invoke this method use method='pape'.
#' 2) We provide an improved compound Poisson approximation that
#' uses more appropriate statistical assumptions concerning
#' overlapping motif hits and that can be used with order-d
#' background models as well. The improved version is used by default
#' with method='kopp'.
#' Note: Only method='kopp' supports the computation
#' of the distribution of the number of motif hits w.r.t. scanning
#' a single DNA strand (see \code{\link{probOverlapHit}}).
#'
#' An Overlap-object is created by 
#' \code{\link{probOverlapHit}}
#' In contrast to \code{\link{combinatorialDist}}, this function
#' supports variable-length DNA sequence.
#'
#' @param seqlen Integer-valued vector that defines the lengths of the
#' individual sequences. For a given DNAStringSet, 
#' this information can be retrieved using \code{\link{numMotifHits}}.
#' @inheritParams overlapValid
#' @param method String that defines which method shall be invoked: 'pape' or
#' 'kopp' (see description). Default: method='kopp'.
#'
#' @return List containing
#' \describe{
#' \item{dist}{Distribution of the number of hits}
#' }
#' @examples
#'
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' alpha=0.001
#' gran=0.1
#' motifcounterOption(alpha, gran)
#'
#' # estimate a background model
#' bg=readBackground(seqs,1)
#'
#' # load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the distribution for scanning the forward DNA strand
#' # of 100 individual 150 bp sequences
#' seqlen=rep(150,100)
#'
#' # Compute overlapping probabilities
#' op=motifcounter:::probOverlapHit(motif,bg,singlestranded=TRUE)
#'
#' # Computes  the distribution of the number of motif hits
#' dist=motifcounter:::compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#'
#' # Proceed similarly for scanning both DNA strands
#' op=motifcounter:::probOverlapHit(motif,bg,singlestranded=FALSE)
#' dist=motifcounter:::compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#'
#' @seealso \code{\link{combinatorialDist}}
#' @seealso \code{\link{probOverlapHit}}
#' @seealso \code{\link{numMotifHits}}
compoundPoissonDist=function(seqlen, overlap,method="kopp") {
    overlapValid(overlap)
    # for all practical purposes, a maximal clump size of 60
    # should be enough
    maxclumpsize=60
    # determine the max number of hits automatically from the
    # given sequence length.
    # even though it is a little bit of computational overhead,
    # set the maxhits to the total sequence length
    maxhits=sum(seqlen)
    dist=numeric(maxhits+1)
    if (method=="kopp") {
        res=.C("motifcounter_compoundPoisson_useBeta", overlap$alpha,
            overlap$beta, overlap$beta3p, overlap$beta5p,
            as.numeric(dist), as.integer(length(seqlen)),
            as.integer(seqlen),
            as.integer(maxhits), as.integer(maxclumpsize),
            length(overlap$beta),
            as.integer(overlap$singlestranded),PACKAGE="motifcounter")
        dist=res[[5]]
    } else if (method=="pape") {
        if (overlap$singlestranded==TRUE) {
            stop("The Pape et al. implementation of the compound
                Poisson distribution can only be for scanning both DNA strands.
                Use probOverlapHit(singlestranded=F).")
        }
        res=.C("motifcounter_compoundPoissonPape_useGamma", overlap$gamma,
            as.numeric(dist), as.integer(length(seqlen)), as.integer(seqlen),
            as.integer(maxhits), as.integer(maxclumpsize),
            length(overlap$beta),
            PACKAGE="motifcounter")
        dist=res[[2]]
    } else {
        stop("The method must be 'kopp' or 'pape'")
    }
    return (list(dist=dist))
}
