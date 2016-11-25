#' Set parameters for the enrichment analysis
#' 
#' This function sets some parameters that are necessary for
#' the computations of the `mdist` package.
#' 
#' 
#' @param alpha Significance level for calling
#' motif hits. E.g. alpha=0.001 amounts to calling 1 in 1000 positions
#' in the DNA motif hit by chance if a single DNA strand is scanned.
#' @param gran The score granularity which is used
#' for discretizing the score range. Decreasing gran values 
#' result an increased number
#' of discrete bins that represent the real-valued score range
#' and consequently mmomre accurate score distribution approximations.
#' On the other hand, decreasing gran incurs higher computations costs.
#' Default: gran=0.1
#' @param ncores Number of cores used for parallel processing if openMP is 
#' available. Default: ncores=1.
#' 
#' @return None
#' @examples
#' 
#' 
#' mdistOption(alpha=0.001)
#' 
#' @export
mdistOption=function(alpha, gran=0.1, ncores=1) {
    dummy=.C("mdist_option", 
        as.numeric(alpha), 
        as.numeric(gran),as.integer(ncores),PACKAGE="mdist")
}
