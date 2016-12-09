#' Set parameters for the enrichment analysis
#'
#' This function sets some parameters that are necessary for
#' the computations of the `motifcounter` package.
#'
#'
#' @param alpha Significance level for calling
#' motif hits. E.g. alpha=0.001 amounts to calling one motif hit per strand
#' by chance in a sequence of length 1000.
#' @param gran The score granularity which is used
#' for discretizing the score range. Decreasing gran values
#' will increase number
#' of discrete bins that represent the real-valued score range.
#' This will yield more a accurate score distribution due to less
#' discretization noise, however, it incurs an increase of the computational
#' runtime.  Default: gran=0.1
#' @param ncores Number of cores used for parallel processing if openMP is
#' available. Default: ncores=1.
#'
#' @return None
#' @examples
#'
#'
#' motifcounterOption(alpha=0.001)
#'
#' @export
motifcounterOption=function(alpha, gran=0.1, ncores=1) {
    if (!is.numeric(alpha)) {
        stop("alpha must be numeric")
    }
    if (alpha<=0.0 || alpha>1.0) {
        stop("alpha must be a probability")
    }
    if (!is.numeric(gran)) {
        stop("gran must be numeric")
    }
    if (gran<=0.0) {
        stop("gran must be positive")
    }
    dummy=.C("motifcounter_option",
        as.numeric(alpha),
        as.numeric(gran),as.integer(ncores),PACKAGE="motifcounter")
}

#' Retrieve the false positive level
#'
#' This function returns the current false positive level for
#' calling motif hits in random sequences.
#'
#'
#'
#' @return False positive motif hit probability
#'
sigLevel=function() {
    alpha=0.0
    ret=.C("motifcounter_siglevel",alpha)
    return(ret[[1]])
}
