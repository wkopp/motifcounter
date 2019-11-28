#' Overlap class definition
#'
#' Objects of this class serve as a container that holds
#' parameters for the overlapping hit probabilities.
#'
#' An Overlap object is constructed via the \code{\link{probOverlapHit}}
#'
#' @slot alpha Scalar numeric significance level to call motif matches
#' @slot beta Numeric vector of 
#'              principal overlapping hit probabilities on the same strand.
#' @slot beta3p Numeric vector of
#'              principal overlapping hit probabilities with 3'-overlap.
#' @slot beta5p Numeric vector of
#'              principal overlapping hit probabilities with 5'-overlap.
#' @slot gamma Numeric vector of marginal overlapping hit probabilities.
#' @slot singlestranded logical flag to indicate whether one
#'              or both strands are scanned for motif matches.
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

#' Accessor to slot alpha
#'
#' @param obj An Overlap object
#'
#' @return alpha slot
getAlpha = function(obj) {
    return(obj@alpha)
}

#' Accessor to slot beta
#'
#' @param obj An Overlap object
#'
#' @return beta slot
getBeta = function(obj) {
    return(obj@beta)
}

#' Accessor to slot beta
#'
#' @param obj An Overlap object
#'
#' @return beta5p slot
getBeta5p = function(obj) {
    return(obj@beta5p)
}

#' Accessor to slot beta3p
#'
#' @param obj An Overlap object
#'
#' @return beta3p slot
getBeta3p = function(obj) {
    return(obj@beta3p)
}

#' Accessor to slot gamma
#'
#' @param obj An Overlap object
#'
#' @return gamma slot
getGamma = function(obj) {
    return(obj@gamma)
}

#' Accessor to slot singlestranded
#'
#' @param obj An Overlap object
#'
#' @return singlestranded slot
getSinglestranded = function(obj) {
    return(obj@singlestranded)
}

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
    scorethreshold = scoreThreshold(pfm, bg)

    if (singlestranded == TRUE) {
        res = .C(
            motifcounter_overlapSingleStranded,
            as.numeric(pfm),
            nrow(pfm),
            ncol(pfm),
            as.numeric(alpha),
            as.numeric(beta),
            as.numeric(beta3p),
            as.numeric(beta5p),
            as.numeric(gamma),
            as.numeric(getStation(bg)),
            as.numeric(getTrans(bg)),
            as.integer(getOrder(bg))
        )
    } else {
        res = .C(
            motifcounter_overlap,
            as.numeric(pfm),
            nrow(pfm),
            ncol(pfm),
            as.numeric(alpha),
            as.numeric(beta),
            as.numeric(beta3p),
            as.numeric(beta5p),
            as.numeric(gamma),
            as.numeric(getStation(bg)),
            as.numeric(getTrans(bg)),
            as.integer(getOrder(bg))
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


