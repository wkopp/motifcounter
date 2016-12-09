#' Score distribution
#'
#' This function computes the score distribution for the provided PFM and
#' background. The Score distribution is computed based on an efficient
#' dynamic programming algorithm.
#'
#' @param pfm A position frequency matrix
#' @param bg A Background object
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
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load the motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' dp=scoreDist(motif,bg)
#'
#' @export
scoreDist=function(pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)

    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.integer(scorerange),
        bg$station,bg$trans,as.integer(bg$order),
        PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    ret=.C("motifcounter_scoredist",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(scores),
        as.numeric(dist),
        bg$station,bg$trans,as.integer(bg$order),
        PACKAGE="motifcounter")
    return(list(score=ret[[4]], probability=ret[[5]]))
}


#' Score distribution
#'
#' This function computes the score distribution for a given PFM and
#' a background model. The result is identical to \code{\link{scoreDist}},
#' however, the method employs a less efficient algorithm that
#' enumerates all DNA sequences of the length of the motif.
#' This function is only used for debugging and testing purposes
#' and might require substantial computational
#' resources for long motifs. Therefore, use \code{\link{scoreDist}} instead.
#'
#' @param pfm A position frequency matrix
#' @param bg A Background object
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
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' dp=scoreDistBf(motif,bg)
#'
#' @export
scoreDistBf=function(pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)

    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.numeric(pfm),nrow(pfm),ncol(pfm),
                as.integer(scorerange),
                bg$station,bg$trans,as.integer(bg$order),
                PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    ret=.C("motifcounter_scoredist_bf",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        as.numeric(scores),
        as.numeric(dist),bg$station,bg$trans,as.integer(bg$order),
        PACKAGE="motifcounter")
    return(list(score=ret[[4]], probability=ret[[5]]))
}

#' Score observations
#'
#' This function computes the per-position and per-strand 
#' score in a given DNA sequence
#'
#' @param pfm A position frequency matrix
#' @param seq DNAString
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{fscores}{Vector of scores on the forward strand}
#' \item{rscores}{Vector of scores on the reverse strand}
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
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the per-position and per-strand scores
#' scoreSequence(seqs[[1]],motif,bg)
#'
#' @export
scoreSequence=function(seq,pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)
    if (class(seq)!="DNAString") {
        stop("seq must be a DNAString object")
    }
    if (length(seq)<ncol(pfm)) {
        stop("length(seq) must be at least as long as the motif")
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
            as.integer(slen),bg$station,bg$trans,as.integer(bg$order),
            PACKAGE="motifcounter")
        return(list(fscores=ret[[5]],rscores=ret[[6]]))
    }
}

#' Score profile across multiple sequences
#'
#' This function computes the per-position and per-strand 
#' average score profiles across a set of DNA sequences.
#'
#' @param pfm A position frequency matrix
#' @param seqs DNAStringSet
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{fscores}{Vector of per-position average forward strand scores}
#' \item{rscores}{Vector of per-position average reverse strand scores}
#' }
#'
#' @examples
#'
#'
#'
#' # Load the DNA sequence and a Motif
#' seqfile=system.file("extdata","oct4_chipseq.fa", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' seqs=seqs[1:10]
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score profile
#' scoreSequenceProfile(seqs,motif,bg)
#'
#' @export
scoreSequenceProfile=function(seqs,pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)
    if (class(seqs)!="DNAStringSet") {
        stop("seq must be a DNAString object")
    }
    if (any(lenSequences(seqs)!=length(seqs[[1]]))) {
        stop("all sequences must be equally long")
    }

    fscores=sapply(seqs, function(seq,pfm,bg) {
        scoreSequence(seq,pfm,bg)$fscores}, 
        pfm,bg)
    fscores=apply(fscores,1,mean)

    rscores=sapply(seqs, function(seq,pfm,bg) {
        scoreSequence(seq,pfm,bg)$rscores}, 
        pfm,bg)
    rscores=apply(rscores,1,mean)
    return (list(fscores=fscores,rscores=rscores))
}

#' Score histogram on a single sequence
#'
#' This function computes the score histogram
#' for a given sequence
#'
#'
#' @param seq DNAString
#' @param pfm A position frequency matrix
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{score}{Vector of score bins}
#' \item{probability}{Empirical frequencies associated with each score}
#' }
#'
#'
scoreHistogramSingleSeq=function(seq,pfm, bg) {
    motifValid(pfm)
    backgroundValid(bg)
    if (class(seq)!="DNAString") {
        stop("seq must be a DNAString or a DNAStringSet object")
    }
    if (length(seq)<ncol(pfm)) {
        stop("length(seq) must be at least as long as the motif")
    }
    scorerange=integer(1)
    scorerange=.C("motifcounter_scorerange",
                as.numeric(pfm),nrow(pfm),ncol(pfm),
                as.integer(scorerange),
                bg$station,bg$trans,as.integer(bg$order),
                PACKAGE="motifcounter")[[4]]
    scores=numeric(scorerange); dist=numeric(scorerange)
    length(scores)

    ret=.C("motifcounter_scorehistogram",
        as.numeric(pfm),nrow(pfm),ncol(pfm),
        toString(seq),as.integer(length(seq)),
        as.numeric(scores), as.numeric(dist),
        bg$station,bg$trans,as.integer(bg$order),
        PACKAGE="motifcounter")
    result=list(score=ret[[6]], frequency=ret[[7]])

    return(result)
}


#' Score histogram
#'
#' This function computes the score histogram
#' for a given sequence
#'
#'
#' @param pfm A position frequency matrix
#' @param seq DNAString or DNAStringSet
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{score}{Vector of score bins}
#' \item{probability}{Empirical frequencies associated with each score}
#' }
#' @examples
#'
#' # Set the the significance level and the score granularity
#' motifcounterOption(alpha=0.01, gran=0.1)
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # generate the simulated score distribution on
#' # sequences of length 1kb using 1000 samples
#' seq=generateDNAString(1000,bg)
#'
#' # Compute the empirical score histogram
#' scoreHistogram(seqs,motif,bg)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreHistogram=function(seq,pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)

    if (class(seq)=="DNAString") {
        result=scoreHistogramSingleSeq(seq,pfm,bg)
    } else if (class(seq)=="DNAStringSet") {
        his=lapply(seq, scoreHistogramSingleSeq,pfm,bg)
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
#' @param bg A Background object
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
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score threshold
#' scoreThreshold(motif,bg)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreThreshold=function(pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)
    scoredist=scoreDist(pfm,bg)

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
