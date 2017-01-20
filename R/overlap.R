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
#' @return A list containing various overlapping hit probabilities.
#' The list contains the following entries
#'     \describe{
#'     \item{alpha}{False positive probability}
#'     \item{beta}{Vector of overlapping hit probability for hits
#'             occurring on the same strand. Each element corresponds to
#'             relative distance of the starts of the two hits.}
#'     \item{beta3p}{Vector of overlapping hit probability for a
#'             forward strand hit that is followed by a reverse strand
#'             hit. Each element corresponds to
#'             relative distance of the starts of the two hits.}
#'     \item{beta5p}{Vector of overlapping hit probability for a
#'             reverse strand hit that is followed by a forward strand
#'             hit. Each element corresponds to
#'             relative distance of the starts of the two hits.}
#'     \item{gamma}{Vector of overlapping hit probabilities across
#'             all configurations. In contrast to beta, beta3p and beta5p,
#'             gamma is not corrected for intermediate motif hit events.}
#'     \item{singlestranded}{singlestranded flag that is prescribed 
#'             as input argument.}
#'     }
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
    backgroundValid(bg)
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
            as.numeric(bg$station),
            as.numeric(bg$trans),
            as.integer(bg$order),
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
            as.numeric(bg$station),
            as.numeric(bg$trans),
            as.integer(bg$order),
            PACKAGE = "motifcounter"
        )
    }
    overlap = list(
        alpha = res[[4]],
        beta = res[[5]],
        beta3p = res[[6]],
        beta5p = res[[7]],
        gamma = res[[8]],
        singlestranded = singlestranded
    )
    class(overlap) = "Overlap"
    return (overlap)
}


#' Check validity of Overlap
#'
#' This function checks if the Overlap object is valid. The function throws
#' an error if the object does not represent an Overlap object.
#'
#' @param overlap An Overlap object which is created by
#' \code{\link{probOverlapHit}}.
#' @return None
#'
#' @examples
#' 
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs,1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute overlapping hit probabilities for scanning both DNA strands
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#' 
#' # Check valididity
#' motifcounter:::overlapValid(op)
#'
overlapValid = function(overlap) {
    if (class(overlap) != "Overlap") {
        stop("overlap must be an Overlap object.
            Use probOverlapHit() to construct one.")
    }
    if (length(overlap) != 6) {
        stop(
            "Overlap object inconsistent.
            overlap must contain 6 elements.
            Use 'probOverlapHit' to construct overlap properly."
        )
    }
}
