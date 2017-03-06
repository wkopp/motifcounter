#' Background class definition
#'
#' Objects of this class serve as a container that holds
#' parameters for the Background model.
#'
#' A Background model is constructed 
#' via \code{\link{readBackground}}.
#'
#' @slot station Stationary probabilities
#' @slot trans Transition probabilities
#' @slot counts k-mer counts
#' @slot order Background model order
.Background <- setClass("Background",
    slots = c(
                station = "numeric",
                trans = "numeric",
                counts = "integer",
                order = "integer"
            ),
    validity = function(object) {
        msg <- character(0)
        # check slot dimensions
        if (length(object@station) != max(4, 4 ^ (object@order))) {
            msg <- paste(strwrap("Inconsistent Background object.
                Use readBackground() to construct a Background object."),
                collapse = "\n")
        }
        if (!length(msg) && length(object@trans) != 4 ^ (object@order + 1)) {
            msg <- paste(strwrap("Inconsistent Background object.
                Use readBackground() to construct a Background object."),
                collapse = "\n")
        }
        # check if slots are normalized
        if (!length(msg) && !isTRUE(all.equal(sum(object@station), 1))) {
            msg <- paste(strwrap("Inconsistent Background object.
                Use readBackground() to construct a Background object."),
                collapse = "\n")
        }
        if (!length(msg) && 
            !isTRUE(all.equal(sum(object@trans), 4^object@order))) {
            msg <- paste(strwrap("Inconsistent Background object.
                Use readBackground() to construct a Background object."),
                collapse = "\n")
        }
        if (length(msg)) msg else TRUE
    }
)



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
    stopifnot (is(seqs, "DNAStringSet"))
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
    background = .Background(
        station = dummy[[2]],
        trans = dummy[[3]],
        counts = as.integer(dummy[[1]]),
        order = as.integer(order)
    )
    return (background)
}


