#' Score distribution
#'
#' This function computes the score distribution for the provided PFM and
#' background. The Score distribution is computed based on an efficient
#' dynamic programming algorithm.
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

#' Score observations
#'
#' This function computes the observed score in a given DNA sequence
#'
#' @param pfm A position frequency matrix
#' @param seq DNAString
#' @return List containing
#' \describe{
#' \item{fscores}{Scores on the forward strand}
#' \item{rscores}{Scores on the reverse strand}
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
#' seq=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load the order-1 background model from the DNA sequence
#' readBackground(seqfile,1)
#'
#' # Load the motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' scoreSequence(motif,seq[[1]])
#'
#' @export
scoreSequence=function(pfm,seq) {
    motifValid(pfm)
    if (class(seq)!="DNAString") {
        stop("seq must be a DNAString object")
    }
    # allocate the scores per position
    slen=length(seq)-ncol(pfm)+1

    # Dirty at the moment. Init scores with very small values
    fscores=rep(-1e10,slen)
    rscores=rep(-1e10,slen)

    if (Biostrings::countPattern("N",seq)>0 ||
        Biostrings::countPattern("n",seq)>0) {
        return(list(fscores=fscores,rscores=rscores))
    } else {
        ret=.C("motifcounter_scoresequence",
            as.numeric(pfm),nrow(pfm),ncol(pfm),toString(seq),
            as.numeric(fscores),as.numeric(rscores),
            as.integer(slen),PACKAGE="motifcounter")
        return(list(fscores=ret[[5]],rscores=ret[[6]]))
    }
}


#' Score histogram on a single sequence
#'
#' This function computes the observed score histogram
#' for a given sequence
#'
#'
#' @param seq DNAString
#' @param pfm A position frequency matrix
#' @return List containing
#' \describe{
#' \item{score}{Vector of score bins}
#' \item{probability}{Frequencies}
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
#' readBackgroundForSampling(seqfile,1)
#'
#' # Load a motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' seq=generateDNAString(1000)
#'
#' # generate the simulated score distribution on
#' # sequences of length 1kb using 1000 samples
#' scoreHistogram(motif,seq)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreHistogramSingleSeq=function(seq,pfm) {
    motifValid(pfm)

    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.numeric(pfm),nrow(pfm),ncol(pfm),
                as.integer(scorerange),PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    length(scores)
    if (class(seq)!="DNAString") {
        stop("seq must be a DNAString or a DNAStringSet object")
    }

    ret=.C("motifcounter_scorehistogram",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        toString(seq),as.integer(length(seq)),
        as.numeric(scores), as.numeric(dist),
        PACKAGE="motifcounter")
    result=list(score=ret[[6]], frequency=ret[[7]])

    return(result)
}


#' Score histogram
#'
#' This function computes the observed score histogram
#' for a given sequence
#'
#'
#' @param pfm A position frequency matrix
#' @param seq DNAString or DNAStringSet
#' @return List containing
#' \describe{
#' \item{score}{Vector of score bins}
#' \item{probability}{Frequencies}
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
#' readBackgroundForSampling(seqfile,1)
#'
#' # Load a motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' seq=generateDNAString(1000)
#'
#' # generate the simulated score distribution on
#' # sequences of length 1kb using 1000 samples
#' scoreHistogram(motif,seq)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreHistogram=function(pfm,seq) {
    motifValid(pfm)

    if (class(seq)=="DNAString") {
        result=scoreHistogramSingleSeq(seq,pfm)
    } else if (class(seq)=="DNAStringSet") {
        his=lapply(seq, scoreHistogramSingleSeq,pfm)
        nseq=length(his)
        scores=his[[1]]$score
        nrange=length(his[[1]]$frequency)
        his=lapply(his,function(x) {x$frequency})
        his=unlist(his)
        his=matrix(his,nrange,nseq)
        his=apply(his,1,sum)
        freq=his

        result=list(score=scores, frequency=freq)
    } else {
        stop("seq must be a DNAString or a DNAStringSet object")
    }

    return(result)
}

#' Score threshold
#'
#' This function computes the score threshold for a desired
#' false positive rate to obtain a motif hit.
#' The returned false positive rate usually differs slightly
#' from the one that is set in \code{\link{motifcounterOption}}, because
#' of the discrete nature of the sequences.
#'
#' @param pfm A position frequency matrix
#' @return List containing
#' \describe{
#' \item{threshold}{Score threshold}
#' \item{alpha}{False positive probability}
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
#' scoreThreshold(motif)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreThreshold=function(pfm) {
    motifValid(pfm)

    scoredist=scoreDist(pfm)

    # find quantile
    ind=which(1-cumsum(scoredist$probability)<=sigLevel())
    if (length(ind)<=1) {
        stop("The significance level is too stringent for the given motif.
            There won't be any motif hits at that level.
            Use motifcounterOption() to prescribe a less stringent alpha")
    }
    ind=ind[2:length(ind)]
    alpha=sum(scoredist$probability[ind])

    ind=min(ind)
    threshold=scoredist$score[ind]

    return(list(threshold=threshold, alpha=alpha))
}
