#' Empirical score distribution
#' 
#' This function estimates the empirical score distribution 
#' by simulating random DNA sequences based on the 
#' background model. Subsequently, the random sequences are scanned
#' with the scoring measure which yields the score distribution.
#' 
#' 
#' @param pfm A position frequency matrix
#' @param seqlen Length of the sequence to be
#' simulated.
#' @param nsim Number of samples to be generated.
#' @return List containing 
#' \describe{
#' \item{score}{Vector of scores}
#' \item{probability}{Vector of score probabilities}
#' }
#' @examples
#' 
#' # Set the the significance level and the score granularity
#' motifcounterOption(alpha=0.01, gran=0.1)
#' 
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' 
#' # Load an order-1 background model
#' readBackground(seqfile,1)
#' 
#' # Load a motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' # generate the simulated score distribution on 
#' # sequences of length 1kb using 1000 samples 
#' simulateScoreDist(motif,seqlen= 1000,nsim=1000)
#' 
#' @seealso \code{\link{scoreDist}}
#' @export
simulateScoreDist=function(pfm,seqlen, nsim) {
    motifValid(pfm)
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.numeric(pfm),nrow(pfm),ncol(pfm),
                as.integer(scorerange),PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    length(scores)
    ret=.C("motifcounter_simulateScores", 
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(scores), 
        as.numeric(dist), as.integer(seqlen), 
        as.integer(nsim),PACKAGE="motifcounter")
    return(list(score=ret[[4]], probability=ret[[5]]))
}



#' Empirical number of motif hits distribution
#' 
#' This function simulates random DNA sequences according to 
#' the background model and
#' subsequently counts how many motif hits occur in them.
#' Doing that repeatedly eventually established the empirical
#' distribution of the number of motif hits.
#' 
#' @param pfm A position frequency matrix
#' @param seqlen Integer-valued vector contanting the individual
#' sequence lengths.
#' @param maxhits Maximal number of hits.
#' @param nsim Number of random samples.
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits.
#' @return Empirical number of motif hits distribution
#' @examples
#' 
#' 
#' motifcounterOption(0.01, 0.01)
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' 
#' readBackground(seqfile,1)
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' seqlen=rep(150,100)
#' simc=simulateNumHitsDist(motif, seqlen,
#'         maxhits=1000,nsim=100,singlestranded=FALSE)
#' 
#' simc=simulateNumHitsDist(motif, seqlen,
#'         maxhits=1000,nsim=100,singlestranded=TRUE)
#' 
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
#' @export
simulateNumHitsDist=function(pfm,seqlen, maxhits, nsim, singlestranded=FALSE) {
    motifValid(pfm)
    if (length(seqlen)<=0) {
        stop("seqlen must be non-empty")
    }
    dist=numeric(maxhits+1);
    ret=.C("motifcounter_simulateCountDistribution", 
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(dist), 
        as.integer(nsim), as.integer(length(seqlen)),
        as.integer(seqlen), as.integer(maxhits), 
        as.integer(singlestranded),PACKAGE="motifcounter")
    return(list(dist=ret[[4]]))
}

