#' Motif overlap probabilities
#' 
#' This function computes the self-overlapping probabilites of the motif for
#' the loaded PWM and background based on dynamic programming.
#' 
#' 
#' @param singlestranded Boolian which defines whether the overlapping hit
#' probabilities shall be computed with respect to scanning both DNA strands or
#' only one strand.  Default: Both strands are scanned for DNA sequences.
#' @return A list containing various overlapping hit probabilities,  including

#'             beta ... the overlapping hit probability for hits occurring on the same strand
#'              
#' @examples
#' 
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' alpha=0.001
#' gran=0.1
#' mdistOption(alpha, gran)
#' 
#' # estimate background model from seqfile
#' readBackground(seqfile,1)
#' 
#' # load motif model from pwmfile
#' readMotif(pwmfile,0.01)
#' 
#' # compute the overlap probabilities for scanning both DNA strands
#' op=probOverlapHit(singlestranded=FALSE)
#' 
#' # compute the overlap probabilities for scanning a single DNA strand
#' op=probOverlapHit(singlestranded=TRUE)
#' 
#' @export
probOverlapHit=function(singlestranded=FALSE) {
    mlen=motifLength()
    alpha=numeric(1)
    beta=numeric(mlen)
    beta3p=numeric(mlen)
    beta5p=numeric(mlen)
    gamma=numeric(3*mlen)
    if (singlestranded==TRUE) {
        res=.C("mdist_overlapSingleStranded", alpha, beta, beta3p, beta5p, gamma,PACKAGE="mdist")
    } else {
        res=.C("mdist_overlap", alpha, beta, beta3p, beta5p, gamma,PACKAGE="mdist")
    }
    return (list(alpha=res[[1]],beta=res[[2]],beta3p=res[[3]],beta5p=res[[4]],
  		 gamma=res[[5]], singlestranded=singlestranded))
}



#' Compound Poisson Approximation
#' 
#' This function computes the distribution of the number of motif hits in a set
#' of random DNA sequences according to a compound Poisson approximation.
#' 
#' 
#' @param seqlen Vector of integers that contains the individual sub-sequence
#' lengths. E.g. seqlen=c(100,200) means that the set of sequences is made up
#' of a sequence with 100 bp and 200 bp. Note that for the 'combinatorial
#' model' all sub-sequences must be of equal length.
#' @param overlap Overlapping hit probabilities. This argument is the result of
#' evaluating the method 'probOverlapHit', which is essential for both, the
#' compound Poisson model and the combinatorial model.
#' @param maxhits Integer that sets the maximal number of motif hits. By
#' default, maxhits is set to 1000 motif hits. In case, the user expects to
#' find more than 1000 hits by chance in the sequence, this parameter needs to
#' be increased.
#' @param maxclumpsize Integer that defines the maximal size of a clump.
#' Default: maxclumpsize=60
#' @param method String that defines which method shall be invoked 'pape' or
#' 'kopp'. Default: method='kopp'
#' @examples
#' 
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' alpha=0.001
#' gran=0.1
#' mdistOption(alpha, gran)
#' 
#' # estimate background model from seqfile
#' readBackground(seqfile,1)
#' 
#' # load motif model from pwmfile
#' readMotif(pwmfile, 0.01)
#' 
#' # compute the distribution for scanning a single DNA strand
#' 
#' #Compute overlapping probabilities
#' op=probOverlapHit(singlestranded=TRUE)
#' 
#' # Computes the distribution of the number of motif hits
#' seqlen=rep(150,100)
#' dist=compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#' 
#' # compute the distribution for scanning both DNA strands
#' 
#' #Compute overlapping probabilities
#' op=probOverlapHit(singlestranded=FALSE)
#' 
#' # Computes the distribution of the number of motif hits
#' seqlen=rep(150,100)
#' dist=compoundPoissonDist(seqlen, op)
#' #plot(1:length(dist$dist)-1, dist$dist)
#' 
#' @export
compoundPoissonDist=function(seqlen, overlap, 
  maxhits=1000, maxclumpsize=60, method="kopp") {
  dist=numeric(maxhits+1)
  if (method=="kopp") {
    res=.C("mdist_compoundPoisson_useBeta", overlap$alpha,
        overlap$beta, overlap$beta3p, overlap$beta5p,
        as.numeric(dist), as.integer(length(seqlen)),
        as.integer(seqlen),
        as.integer(maxhits), as.integer(maxclumpsize),as.integer(overlap$singlestranded),PACKAGE="mdist")
    dist=res[[5]]
  } else if (method=="pape") {
    if (overlap$singlestranded==TRUE) {
        stop("The Pape implementation of the compound Poisson distribution works only for scanning both DNA strands.
             Use probOverlapHit(singlestranded=F).")
    }
    res=.C("mdist_compoundPoissonPape_useGamma", overlap$gamma,
        as.numeric(dist), as.integer(length(seqlen)), as.integer(seqlen),
        as.integer(maxhits), as.integer(maxclumpsize),PACKAGE="mdist")
    dist=res[[2]]
  } else {
      stop("The method must be 'kopp' or 'pape'")
  }
  return (list(dist=dist))
}

