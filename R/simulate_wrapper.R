#' Empirical score distribution
#'
#' This function estimates the empirical score distribution
#' by simulating random DNA sequences based on the
#' background model. Subsequently, the random sequences are scanned
#' with the scoring measure which yields the score distribution.
#'
#'
#' @param pfm A position frequency matrix
#' @param bg A Background object
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
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # generate the simulated score distribution on
#' # sequences of length 1kb using 1000 samples
#' scoreDistEmpirical(motif,bg,seqlen= 1000,nsim=1000)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreDistEmpirical=function(pfm,bg,seqlen, nsim) {
    motifValid(pfm)
    backgroundValid(bg)
    seqs=generateDNAStringSet(rep(seqlen,nsim),bg)
    sh=scoreHistogram(seqs,pfm,bg)

    probs=sh$frequency/sum(sh$frequency)
    return(list(score=sh$score, probability=probs))
}

#' Generate DNAString
#'
#' This function generates a DNAString
#' by sampling from the the background model.
#'
#' @param len Length of the sequence
#' @param bg A Background object
#'
#' @return DNAString
#
#' @examples
#'
#' # Set the the significance level and the score granularity
#' motifcounterOption(alpha=0.01, gran=0.1)
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # generates a sequences of length 1kb
#' generateDNAString(100,bg)
#'
#' @export
generateDNAString=function(len,bg) {
    len=as.integer(len)
    backgroundValid(bg)
    if (is.na(len)) {
        stop("len must be a positive integer")
    }
    if (len<=0) {
        stop("len must be a positive integer")
    }
    if (len<bg$order) {
        stop(sprintf("len must be at least %d",bg$order))
    }
    seq=paste(rep("N",len),collapse="",sep="")
    ret=.C("motifcounter_generateRndSeq", as.character(seq),
        as.integer(len),
        bg$station,bg$trans,as.integer(bg$order),
        PACKAGE="motifcounter")
    return(Biostrings::DNAString(ret[[1]]))
}

#' Generate DNAStringSet
#'
#' This function generates a DNAStringSet
#' by sampling from the background model.
#'
#'
#' @param len Vector of individual sequence lengths
#' @param bg A Background object
#'
#' @return DNAStringSet
#
#' @examples
#'
#' # Set the the significance level and the score granularity
#' motifcounterOption(alpha=0.01, gran=0.1)
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # generates sequences of various lengths
#' generateDNAStringSet(10:50,bg)
#'
#' @export
generateDNAStringSet=function(len,bg) {
    backgroundValid(bg)
    seqs=c()
    for (i in 1:length(len)) {
        seqs=c(seqs,generateDNAString(len[i],bg))
    }

    return(Biostrings::DNAStringSet(seqs))
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
#' @param bg A Background object
#' @param seqlen Integer-valued vector contanting the individual
#' sequence lengths.
#' @param nsim Number of random samples.
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits.
#' @return Empirical number of motif hits distribution
#' @examples
#'
#'
#' motifcounterOption(0.01, 0.01)
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#'
#' bg=readBackground(seqs,1)
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' seqlen=rep(150,1)
#' simc=simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=FALSE)
#'
#' simc=simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=TRUE)
#'
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
#' @export
simulateNumHitsDist=function(pfm,bg,seqlen, nsim, singlestranded=FALSE) {
    motifValid(pfm)
    backgroundValid(bg)
    if (length(nsim)<=0) {
        stop("nsim must be strictly positive")
    }

    freq=rep(0,100)
    for (i in 1:nsim) {
        seqs=generateDNAStringSet(seqlen,bg)
        nom=numMotifHits(seqs,pfm,bg,singlestranded)
        nom=sum(nom$numofhits)
        if (nom>length(freq)+1) {
            #expand freq
            new_freq=rep(0,nom+1)
            new_freq[1:length(freq)]=freq
            freq=new_freq
        }
        freq[nom+1]=freq[nom+1]+1
    }
    freq=freq/sum(freq)
    return(list(dist=freq))
}
