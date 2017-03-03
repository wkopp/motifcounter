#' Enrichment of motif hits
#'
#' This function determines whether a given motif is enriched in a given
#' DNA sequences. 
#' 
#' Enrichment is tested by comparing the observed 
#' number of motif hits against a theoretical distribution of the number
#' of motif hits in random DNA sequences.
#' Optionally, the theoretical distribution of the number of motif
#' hits can be evaluated by either a 'compound Poisson model' 
#' or the 'combinatorial model'.
#' Additionally, the enrichment test can be conducted with respect
#' to scanning only the forward strand or both strands of the DNA
#' sequences. The latter option is only available for the
#' 'compound Poisson model'
#' @include count_wrapper.R
#'
#' @inheritParams numMotifHits
#' @param method String that defines whether to use
#' the 'compound' Poisson approximation' or the 'combinatorial' model.
#' Default: method='compound'.
#'
#' @return List that contains
#' \describe{
#' \item{pvalue}{P-value for the enrichment test}
#' \item{fold}{Fold-enrichment with respect to the expected number of hits}
#' }
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#' 
#' # Load background
#' bg = readBackground(seqs, 1)
#' 
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # 1 ) Motif enrichment test w.r.t. scanning a *single* DNA strand
#' # based on the 'Compound Poisson model'
#' 
#' result = motifEnrichment(seqs, motif, bg,
#'             singlestranded = TRUE, method = "compound")
#'
#' # 2 ) Motif enrichment test w.r.t. scanning *both* DNA strand
#' # based on the 'Compound Poisson model'
#' 
#' result = motifEnrichment(seqs, motif, bg, method = "compound")
#'
#' # 3 ) Motif enrichment test w.r.t. scanning *both* DNA strand
#' # based on the *combinatorial model*
#'
#' result = motifEnrichment(seqs, motif, bg, singlestranded = FALSE,
#'             method = "combinatorial")
#'
#' @seealso \code{\link{compoundPoissonDist}}, \code{\link{combinatorialDist}}
#' @export
motifEnrichment = function(seqs, pfm, bg, singlestranded = FALSE, 
                            method = "compound") {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    #compute overlapping hit probs
    overlap = probOverlapHit(pfm, bg, singlestranded)

    # detemine the number of motif hits
    observations = numMotifHits(seqs, pfm, bg, singlestranded)
    if (method == "compound") {
        dist = compoundPoissonDist(observations$lseq, overlap)
    } else if (method == "combinatorial") {
        dist = combinatorialDist(observations$lseq, overlap)
    } else {
        stop("Invalid method: 'method' must be 'compound' or 'combinatorial'")
    }
    p = sum(dist$dist[(sum(observations$numofhits) + 1):length(dist$dist)])

    return (list(
        pvalue = p,
        fold = sum(observations$numofhits) /
            sum(dist$dist * seq(0, length(dist$dist) - 1))
    ))
}
