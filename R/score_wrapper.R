#' Score distribution
#' 
#' This function computes the score distribution for the provided PFM and 
#' background. The Score distribution is computed based on dynamic
#' programming.
#' 
#' @param pfm A position frequency matrix
#' @return List containing 
#' \describe{
#' \item{score}{Vector of motif scores}
#' \item{probability}{Score probabilities}
#' }
#' 
#' @examples
#' 
#' 
#' # Set the significance level and granularity for the score computation
#' motifcounterOption(alpha=0.01,gran=0.1)
#' 
#' # Load the DNA sequence and a Motif
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' 
#' # Load the order-1 background model from the DNA sequence
#' readBackground(seqfile,1)
#' 
#' # Load the motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' # Compute the score distribution
#' dp=scoreDist(motif)
#' 
#' @export
scoreDist=function(pfm) {
    motifValid(pfm)
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.integer(scorerange),PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    ret=.C("motifcounter_scoredist",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(scores),
        as.numeric(dist),PACKAGE="motifcounter")
    return(list(score=ret[[4]], probability=ret[[5]]))
}

#' Score distribution
#' 
#' This function computes the score distribution for the provided PFM and 
#' background model. The result is identical to \code{\link{scoreDist}},
#' however, the method employs a less efficient algorithm that
#' enumerates all DNA sequences of the length of the motif. 
#' This function is only used for debugging purposes
#' and might require substantial computational
#' resources for long motifs. Therefore, use \code{\link{scoreDist}} instead.
#' 
#' @param pfm A position frequency matrix
#' @return List containing
#' \describe{
#' \item{score}{Vector of motif scores}
#' \item{probability}{Score probabilities}
#' }
#' 
#' @seealso \code{\link{scoreDist}}
#' @examples
#' 
#' 
#' # Set the significance level and granularity for the score computation
#' motifcounterOption(alpha=0.01,gran=0.1)
#' 
#' # Load the DNA sequence and a Motif
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' 
#' # Load the order-1 background model from the DNA sequence
#' readBackground(seqfile,1)
#' 
#' # Load the motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' # Compute the score distribution
#' dp=scoreDistBf(motif)
#' 
#' @export
scoreDistBf=function(pfm) {
    motifValid(pfm)
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.numeric(pfm),nrow(pfm),ncol(pfm),
                as.integer(scorerange),PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    ret=.C("motifcounter_scoredist_bf",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(scores),
        as.numeric(dist),PACKAGE="motifcounter")
    return(list(score=ret[[4]], probability=ret[[5]]))
}

