#' Compound Poisson Approximation
#'
#' This function approximates the distribution of the number of motif hits
#' that emerges from a random DNA sequence of a given length.
#' 
#' The distribution can be determined in two alternative ways:
#' \enumerate{
#' \item A re-implemented version of the algorithm that was
#' described in Pape et al. \emph{Compound poisson approximation
#' of the number of occurrences of a position
#' frequency matrix (PFM) on both strands.} 2008
#' can be invoked using method='pape'.
#' The main purpose of this implementation concerns 
#' benchmarking an improved approximation.
#' In contrast to the original model, this implementation 
#' can be used with general order-d Markov models.
#' \item We provide an improved compound Poisson approximation that
#' uses more appropriate statistical assumptions concerning
#' overlapping motif hits and that can be used with order-d
#' background models as well. The improved version is used by default
#' with method='kopp'.
#' Note: Only method='kopp' supports the computation
#' of the distribution of the number of motif hits w.r.t. scanning
#' a single DNA strand (see \code{\link{probOverlapHit}}).
#' }
#'
#'
#' @param seqlen Integer-valued vector that defines the lengths of the
#' individual sequences. For a given DNAStringSet, 
#' this information can be retrieved using \code{\link{numMotifHits}}.
#' @inheritParams overlapValid
#' @param method String that defines which method shall be invoked: 'pape' or
#' 'kopp' (see description). Default: method = 'kopp'.
#'
#' @return List containing
#' \describe{
#' \item{dist}{Distribution of the number of hits}
#' }
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
#' # Use 100 individual sequences of length 150 bp each
#' seqlen = rep(150, 100)
#'
#' # Compute overlapping probabilities
#' # for scanning the forward DNA strand only
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = TRUE)
#'
#' # Computes  the compound Poisson distribution
#' dist = motifcounter:::compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#'
#' # Compute overlapping probabilities
#' # for scanning the forward DNA strand only
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#' 
#' # Computes  the compound Poisson distribution
#' dist = motifcounter:::compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#'
#' @seealso \code{\link{combinatorialDist}}
#' @seealso \code{\link{probOverlapHit}}
#' @seealso \code{\link{numMotifHits}}
compoundPoissonDist = function(seqlen, overlap, method = "kopp") {
    overlapValid(overlap)
    
    # Length must be at least as long as the motif
    sl = sum(vapply(seqlen, function(sl, ml) {
        sl - ml + 1
    }, FUN.VALUE = 0,
    ml = length(overlap$beta)))

    if (sl <= 0) {
        return (list(dist = 1))
    }
    
    # for all practical purposes, a maximal clump size of 60
    # should be enough
    maxclumpsize = 60
    
    # determine the max number of hits automatically from the
    # given sequence length.
    # even though it is a little bit of computational overhead,
    # set the maxhits to the total sequence length
    maxhits = sum(seqlen)
    dist = numeric(maxhits + 1)
    if (method == "kopp") {
        res = .C(
            "motifcounter_compoundPoisson_useBeta",
            overlap$alpha,
            overlap$beta,
            overlap$beta3p,
            overlap$beta5p,
            as.numeric(dist),
            as.integer(length(seqlen)),
            as.integer(seqlen),
            as.integer(maxhits),
            as.integer(maxclumpsize),
            length(overlap$beta),
            as.integer(overlap$singlestranded),
            PACKAGE = "motifcounter"
        )
        dist = res[[5]]
    } else if (method == "pape") {
        if (overlap$singlestranded == TRUE) {
            stop(
                paste(strwrap(
                "method = 'pape' does not support single stranded scanning
                for motif hits.
                Please use method = 'kopp' instead."), collapse = "\n")
            )
        }
        res = .C(
            "motifcounter_compoundPoissonPape_useGamma",
            overlap$gamma,
            as.numeric(dist),
            as.integer(length(seqlen)),
            as.integer(seqlen),
            as.integer(maxhits),
            as.integer(maxclumpsize),
            length(overlap$beta),
            PACKAGE = "motifcounter"
        )
        dist = res[[2]]
        } else {
            stop(paste(strwrap("Invalid 'method': The 'method' must 
                 be 'kopp' or 'pape'"), collapse = "\n"))
        }
    return (list(dist = dist))
}
