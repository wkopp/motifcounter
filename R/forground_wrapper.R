#' Check valididity of PFM
#'
#' This function checks if the PFM is valid. The function throws
#' an error if the R matrix does not represent a PFM.
#'
#' @param pfm An R matrix that represents a position frequency matrix
#' @return None
#'
#' @examples
#'
#' # Load motif
#' motiffile = system.file("extdata", "x1.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Check validity
#' motifcounter:::motifValid(motif)
#'
motifValid = function(pfm) {
    stopifnot(is.matrix(pfm))

    #check if matrix has four rows
    stopifnot(nrow(pfm) == 4)

    #check if all entries are positive
    if (!all(pfm > 0)) {
        stop("pfm must be strictly positive. Use 'normalizeMotif'.")
    }
    #check if all columns sum to one
    if (!all(abs(1 - colSums(pfm)) < 0.0000001)) {
        stop("Columns must sum to one. Use 'normalizeMotif'.")
    }
}

#' Normalizes a PFM
#'
#' This function normalizes a PFM and optionally
#' adds pseudo-evidence to each entry of the matrix.
#'
#' @inheritParams motifValid
#' @param pseudo Small numeric pseudo-value that is added
#' to each entry in the PFM in order to ensure strictly positive entries.
#' Default: pseudo = 0.01
#' @return A normalized PFM
#'
#' @examples
#'
#' # Load motif
#' motiffile = system.file("extdata", "x1.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Normalize motif
#' new_motif = normalizeMotif(motif)
#'
#' @export
normalizeMotif = function(pfm, pseudo = 0.01) {
    pfm = pfm + pseudo
    pfm = t(t(pfm) / colSums(pfm))
    return(pfm)
}

#' Check valididity of PFM with background
#'
#' This function checks if the PFM x background combination is valid.
#' The function throws an error if this is not the case.
#'
#' @inheritParams motifValid
#' @template templates
#' @return None
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
#' motiffile = system.file("extdata", "x1.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Check validity
#' motifcounter:::motifAndBackgroundValid(motif, bg)
#'
motifAndBackgroundValid = function(pfm, bg) {
    if (ncol(pfm) < getOrder(bg)) {
        stop(paste(strwrap("ncol(pfm) must be greater or equal
            to the background model order (see order-slot of 
            the Background object)"),
            collapse = "\n"))
    }
}

#' Reverse complements a PFM
#'
#' This function computes the reverse complement of a given PFM.
#'
#' @inheritParams motifValid
#' @return Reverse complemented PFM
#'
#' @examples
#'
#' # Load motif
#' motiffile = system.file("extdata", "x1.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Reverse complement motif
#' revcompmotif = motifcounter:::revcompMotif(motif)
#'
revcompMotif = function(pfm) {
    nrows = nrow(pfm)
    ncols = ncol(pfm)
    cpfm = matrix(rev(as.vector(pfm)), nrow = nrows, ncol = ncols)
    return(cpfm)
}

#' Create position weight representation
#'
#' This function computes the position weight representation for a given
#' PFM and a background model which corresponds to the log-likelihood
#' ratio for each position. For order=0, this corresponds to the
#' commonly used (4 x M) position weight matrix where M denotes
#' the motif length.
#' For higher order background, we obtain a pow(4, order+1) x (M - order + 1)
#' dimensional matrix.
#'
#' @inheritParams motifAndBackgroundValid
#' @return Position weight representation
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
#' motiffile = system.file("extdata", "x1.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Reverse complement motif
#' pwm = motifcounter:::getPositionWeights(motif, bg)
#'
getPositionWeights = function(pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    pwm = .Call(
        motifcounter_getpositionweights,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        getStation(bg),
        getTrans(bg),
        as.integer(getOrder(bg))
    )
    return(matrix(pwm,
                nrow=4^(getOrder(bg) + 1),
                ncol=ncol(pfm) - getOrder(bg)))
}
