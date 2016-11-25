#' Score distribution
#' 
#' This function computes the score distribution for the provided PFM and 
#' background. The Score distribution is computed based on dynamic
#' programming.
#' 
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
#' readMotif(motiffile)
#' 
#' # Compute the score distribution
#' dp=scoreDist()
#' 
#' @export
scoreDist=function() {
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.integer(scorerange),PACKAGE="motifcounter")[[1]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    ret=.C("motifcounter_scoredist",as.numeric(scores),
        as.numeric(dist),PACKAGE="motifcounter")
    return(list(score=ret[[1]], probability=ret[[2]]))
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
#' readMotif(motiffile)
#' 
#' # Compute the score distribution
#' dp=scoreDistBf()
#' 
#' @export
scoreDistBf=function() {
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.integer(scorerange),PACKAGE="motifcounter")[[1]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    .C("motifcounter_scoredist_bf",as.numeric(scores),
        as.numeric(dist),PACKAGE="motifcounter")
}

