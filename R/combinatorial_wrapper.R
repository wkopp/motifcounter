#' Combinatrial distribution for the number of motif hits
#' 
#' This function approxmiates the distribution of the number of motif hits by
#' summing over all combinations of obtaining k hits in a given sequence. The
#' summation is done by employing dynamic programming.
#' 
#' 
#' @param seqlen Is an integer-valued vector that defines the lengths of the
#' individual sequences.
#' @param overlap A list that contains the overlapping hit probabilities which
#' is produced by `probOverlapHit()`
#' @param maxhits Integer that sets the maximal number of motif hits.
#' @return The distribution of the number of motif hits
#' @examples
#' 
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' alpha=0.001
#' mdistOption(alpha)
#' 
#' # estimate background model from seqfile
#' readBackground(seqfile,1)
#' 
#' # load motif model from pwmfile
#' readMotif(pwmfile, 0.01)
#' 
#' # Compute overlap probabilies
#' op=probOverlapHit(singlestranded=FALSE)
#' seqlen=rep(100,1)
#' # Computes the distribution of the number of motif hits
#' dist=combinatorialDist(seqlen, op)
#' 
#' @export
combinatorialDist=function(seqlen, overlap, maxhits=100) {
  if (!all(seqlen==seqlen[1])) {
     stop("The sequences must be of equal length for dynamic programming!");
  }
  if (overlap$singlestranded==TRUE) {
      stop("Currently the combinatorial model only supports scanning
            of both DNA strands")
  }
  dist=numeric(maxhits+1)
  ret=.C("mdist_combinatorialDist", 
        overlap$alpha, overlap$beta, overlap$beta3p, overlap$beta5p,
        as.numeric(dist), as.integer(seqlen[1]),
        as.integer(maxhits), as.integer(length(seqlen)),
        as.integer(overlap$singlestranded),PACKAGE="mdist")
  return(list(dist=ret[[5]]))
}

