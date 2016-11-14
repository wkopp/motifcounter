#' Score distribution
#' 
#' This function computes the score distribution for the provided PFM and order
#' m background. The Score distribution is computed based on dynamic
#' programming.
#' 
#' 
#' @examples
#' 
#' 
#' library(mdist)
#' # Set the significance level and granularity for the score computation
#' mdistOption(alpha=0.01,gran=0.1)
#' 
#' # Load the DNA sequence and a Motif
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab",package="mdist")
#' 
#' # Load the order-1 background model from the DNA sequence
#' readBackground(seqfile,1)
#' 
#' # Load the motif from the pwmfile
#' readMotif(pwmfile)
#' 
#' # Compute the score distribution
#' dp=scoreDist()
#' 
#' @export
scoreDist=function() {
  scorerange=integer(1)
  scorerange=.C("mdist_scorerange",as.integer(scorerange),PACKAGE="mdist")[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  .C("mdist_scoredist",as.numeric(scores),as.numeric(dist),PACKAGE="mdist")
}

#' Score distribution
#' 
#' This function computes the score distribution for the provided PFM and order
#' m background based on enumerating all DNA sequences of the 
#' length of the motif. This function is only used for debugging
#' \code{\link{scoreDist}}, and might require substantial computational
#' resources for long motifs.
#' 
#' 
#' @examples
#' 
#' 
#' library(mdist)
#' # Set the significance level and granularity for the score computation
#' mdistOption(alpha=0.01,gran=0.1)
#' 
#' # Load the DNA sequence and a Motif
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab",package="mdist")
#' 
#' # Load the order-1 background model from the DNA sequence
#' readBackground(seqfile,1)
#' 
#' # Load the motif from the pwmfile
#' readMotif(pwmfile)
#' 
#' # Compute the score distribution
#' dp=scoreDistBf()
#' 
#' @export
scoreDistBf=function() {
  scorerange=integer(1)
  scorerange=.C("mdist_scorerange",as.integer(scorerange),PACKAGE="mdist")[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  .C("mdist_scoredist_bf",as.numeric(scores),as.numeric(dist),PACKAGE="mdist")
}

