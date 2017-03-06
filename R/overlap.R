#' Overlap class definition
#'
#' Objects of this class serve as a container that holds
#' parameters for the Overlap hit probabilities, including the
#' scalar significance level alpha, vectors of principal
#' overlapping hit probabilities (beta, beta3p and beta5p),
#' a vector of marginal overlapping hit probabilities (gamma)
#' and a logical-valued singlestranded flag to decide whether one
#' or both strands are scanned for motif hits.
#' An Overlap object is constructed via the probOverlapHit-constructor function 
#' (see below).
.Overlap <- setClass("Overlap",
    slots = c(
                alpha = "numeric",
                beta = "numeric",
                beta3p = "numeric",
                beta5p = "numeric",
                gamma = "numeric",
                singlestranded = "logical"
            )
)

#' Overlapping motif hit probabilities
#'
#' This function computes a set of self-overlapping probabilites for a
#' motif and background model.
#' 
#' The `gamma`s are determined based on two-dimensional score 
#' distributions (similar as described in Pape et al. 2008),
#' however, they are computed based on an order-d background model.
#' On the other hand, the `beta`s represent overlapping hit probabilities
#' that were corrected for intermediate hits.
#' 
#' @inheritParams numMotifHits
#' @return An Overlap object
#'
#' @examples
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
#' # Compute overlapping hit probabilities for scanning both DNA strands
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#'
#' # Compute overlapping hit probabilities for scanning a single DNA strand
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = TRUE)
#'
probOverlapHit = function(pfm, bg, singlestranded = FALSE) {
    #check if pfm is a matrix
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    alpha = numeric(1)
    beta = numeric(ncol(pfm))
    beta3p = numeric(ncol(pfm))
    beta5p = numeric(ncol(pfm))
    gamma = numeric(3 * ncol(pfm))
    if (singlestranded == TRUE) {
        res = .C(
            "motifcounter_overlapSingleStranded",
            as.numeric(pfm),
            nrow(pfm),
            ncol(pfm),
            as.numeric(alpha),
            as.numeric(beta),
            as.numeric(beta3p),
            as.numeric(beta5p),
            as.numeric(gamma),
            as.numeric(bg@station),
            as.numeric(bg@trans),
            as.integer(bg@order),
            PACKAGE = "motifcounter"
        )
    } else {
        res = .C(
            "motifcounter_overlap",
            as.numeric(pfm),
            nrow(pfm),
            ncol(pfm),
            as.numeric(alpha),
            as.numeric(beta),
            as.numeric(beta3p),
            as.numeric(beta5p),
            as.numeric(gamma),
            as.numeric(bg@station),
            as.numeric(bg@trans),
            as.integer(bg@order),
            PACKAGE = "motifcounter"
        )
    }
    overlap = .Overlap(
        alpha = res[[4]],
        beta = res[[5]],
        beta3p = res[[6]],
        beta5p = res[[7]],
        gamma = res[[8]],
        singlestranded = singlestranded
    )
    return (overlap)
}


