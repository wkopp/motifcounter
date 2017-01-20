#' Estimates a background model from a set of DNA sequences
#'
#' Given a set of DNA sequences and an order, this function
#' estimates an order-d Markov model which is used to characterize
#' random DNA sequences.
#'
#'
#' @inheritParams lenSequences
#' @param order Order of the Markov models that shall be used as the
#' background model. Default: order = 1.
#'
#' @return A Background object
#'
#' @examples
#'
#' # Load sequences
#' file = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(file)
#'
#' # Estimate an order-1 Markov model
#' bg = readBackground(seqs, 1)
#'
#' @export
readBackground = function(seqs, order = 1) {
    stopifnot (class(seqs) == "DNAStringSet")
    stopifnot (order >= 0)
    
    trans = numeric(4 ^ (order + 1))
    
    # collect k-mer frequencies from each individual sequence
    counts = vapply(seqs, function(seq, order, trans) {
        return(
            .C(
                "motifcounter_countfreq",
                toString(seq),
                length(seq),
                trans,
                as.integer(order),
                PACKAGE = "motifcounter"
            )[[3]]
        )
    }, FUN.VALUE = trans, order, trans)
    counts = rowSums(counts)
    if (order == 0) {
        station = numeric(4)
    } else {
        station = numeric(4 ^ order)
    }
    dummy = .C(
        "motifcounter_bgfromfreq",
        as.numeric(counts),
        as.numeric(station),
        as.numeric(trans),
        as.integer(order),
        PACKAGE = "motifcounter"
    )
    background = list(
        station = dummy[[2]],
        trans = dummy[[3]],
        counts = dummy[[1]],
        order = order
    )
    class(background) = "Background"
    return (background)
}

#' Check valididity of Background model
#'
#' This function checks if the Background model is valid. The function throws
#' an error if the object does not represent a Background model.
#'
#' @param bg A Background object
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
#' # Validate Background object
#' motifcounter:::backgroundValid(bg)
#'
backgroundValid = function(bg) {
    if (class(bg) != "Background") {
        stop("bg must be a Background object.
            Use readBackground() to construct bg.")
    }
    if (length(bg) != 4) {
        stop("bg must contain 4 elements.
            Use readBackground() to construct bg.")
    }
    if (length(bg$trans) != 4 ^ (bg$order + 1)) {
        stop("Inconsistent Background object.
            Use readBackground() to construct bg.")
    }
}
