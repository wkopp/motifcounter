#' Set parameters for the enrichment analysis
#'
#' This function sets some global parameters for
#' the `motifcounter` package.
#'
#' alpha=0.001 amounts to calling 
#' one motif hit per strand by chance in a sequence of length 1000 bp. 
#' Decreasing gran will increase number
#' of discrete bins that represent the real-valued score range.
#' This will yield more a accurate score distribution due to less
#' discretization noise, however, it incurs an 
#' increase of the computational burden.
#' 
#' @param alpha Numeric False positive probabililty for calling
#' motif hits by chance. Default: alpha = 0.001
#' @param gran Numeric score granularity which is used
#' for discretizing the score range. Default: gran = 0.1
#' @param ncores Interger number of cores used for parallel processing,
#' if openMP is available. Default: ncores = 1
#'
#' @return None
#' @examples
#'
#' # Prescribe motifcounter Options
#' motifcounterOptions(alpha = 0.001, gran = 0.1, ncores = 1)
#'
#' @export
motifcounterOptions = function(alpha = 0.001, gran = 0.1, ncores = 1) {
    if (alpha <= 0.0 || alpha > 1.0) {
        stop("alpha must be a probability")
    }
    if (gran <= 0.0) {
        stop("gran must be positive")
    }
    if (alpha > 0.05) {
        warning(paste(strwrap(
            "alpha is too chosen high. This might cause
            biases which could lead to false conclusions 
            when using the 'motifEnrichment'.
            Consider prescribing a smaller alpha with 
            'motifcounterOptions()'."),
            collapse = "\n"))
    }
    dummy = .C(
        motifcounter_option,
        as.numeric(alpha),
        as.numeric(gran),
        as.integer(ncores)
    )
}

#' Retrieve the false positive probability
#'
#' This function returns the current false positive level for
#' calling motif hits in random sequences.
#' 
#' The returned value is usually slightly smaller than
#' the prescribed `alpha` in `motifcounterOptions`, because
#' of the discrete nature of sequences.
#'
#'
#' @return False positive probability
#'
#' @examples
#' motifcounter:::sigLevel()
#'
sigLevel = function() {
    alpha = 0.0
    ret = .C(motifcounter_siglevel, alpha)
    return(ret[[1]])
}
