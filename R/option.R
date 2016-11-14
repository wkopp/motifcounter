#' Set parameters for the enrichment analysis
#' 
#' This function sets parameters for the enrichment analysis.
#' 
#' 
#' @param alpha Numerical value that defines the significance level for calling
#' motif hits. E.g. alpha=0.001 amounts to calling 1 in 1000 nucleotides a
#' motif hit by chance if a single DNA strand is scanned.
#' @param gran Numerical value that defines the granularity which shall be used
#' for discretizing the score range. Smaller values of gran might increase the
#' accuracy of the distribution at the cost of an increased computational
#' runtime. Default: gran=0.1
#' @param ncores Number of cores used for parallel processing. Default:
#' ncores=1.
#' @examples
#' 
#' 
#' library(mdist)
#' mdistOption(alpha=0.001)
#' 
#' @export
mdistOption=function(alpha, gran=0.1, ncores=1) {
  dummy=.C("mdist_option", 
           as.numeric(alpha), 
           as.numeric(gran),as.integer(ncores),PACKAGE="mdist")
}
