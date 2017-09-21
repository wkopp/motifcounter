#' Markov model for generating Y_1Y_2_Y3 ...
#'
#' This function implements the Markov model for producing
#' motif matches. The function takes a state probability vector
#' and uses the transition probabilities in order to obtain
#' the state probability at the next time point.
#' This function is used used to determine the stationary distribution
#' of the states.
#'
#' The R interface is only used for the purpose of testing
#' the correctness of the model.
#'
#' @include comppoiss_wrapper.R
#'
#' @inheritParams compoundPoissonDist
#' @param nsteps Number of state transitions to perform
#'
#' @return List containing
#' \describe{
#' \item{dist}{State probability distribution after the given number of steps}
#' }
#' @seealso \code{\link{compoundPoissonDist}}
#' @seealso \code{\link{numMotifHits}}
#' @seealso \code{\link{probOverlapHit}}
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
#' # Compute overlap probabilities
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#'
#'
#' # Computes the state probabilities of the Markov model
#' # (default: after one step)
#' dist = motifcounter:::markovModel(op)
#'
markovModel = function(overlap, nsteps = 1) {
    stopifnot(is(overlap, "Overlap"))

    if (getSinglestranded(overlap) == TRUE) {
        dist = numeric(length(getBeta(overlap)))
        # This might change in the future though
        ret = .C(
            motifcounter_markovmodel_ss,
            as.numeric(dist),
            getAlpha(overlap),
            getBeta(overlap),
            as.integer(nsteps),
            length(getBeta(overlap))
        )
    } else {
        dist = numeric(length(getBeta(overlap))*2 + 2)
        ret = .C(
            motifcounter_markovmodel_ds,
            as.numeric(dist),
            getAlpha(overlap),
            getBeta(overlap),
            getBeta3p(overlap),
            getBeta5p(overlap),
            as.integer(nsteps),
            length(getBeta(overlap))
        )
    }


    return(list(dist = ret[[1]]))
}


#' Computes the Clump start probability based on a Markov model
#'
#' This function leverages a Markov model in order to
#' determine the clump start probability.
#' The computation depends on the selected false positive probability 
#' for calling motif matches 'alpha' and the pre-determined
#' overlapping match probabilities 'beta'.
#'
#' The general idea of the method relies on the fact that
#' for the stationary distribution of the Markov model,
#' motif matches must be observed with probability 'alpha'.
#' Hence, the clump start probability 'tau' is optimized
#' to achieve that goal.
#'
#' The R interface is only used for the purpose of testing
#' the correctness of the model.
#'
#' @include comppoiss_wrapper.R
#'
#' @inheritParams compoundPoissonDist
#'
#' @return Clump start probability 'tau'
#' @seealso \code{\link{compoundPoissonDist}}
#' @seealso \code{\link{numMotifHits}}
#' @seealso \code{\link{probOverlapHit}}
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
#' # Compute overlap probabilities
#' op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)
#'
#'
#' # Computes the clump start probability
#' dist = motifcounter:::computeClumpStartProb(op)
#'
computeClumpStartProb = function(overlap) {
    stopifnot(is(overlap, "Overlap"))

    if (getSinglestranded(overlap) == TRUE) {
        # This might change in the future though
        tau = .Call(
            motifcounter_mcss_check_optimal,
            getAlpha(overlap),
            getBeta(overlap),
            length(getBeta(overlap))
        )
    } else {
        tau = .Call(
            motifcounter_mcds_check_optimal,
            getAlpha(overlap),
            getBeta(overlap),
            getBeta3p(overlap),
            getBeta5p(overlap),
            length(getBeta(overlap))
        )
    }

    return(tau)
}
