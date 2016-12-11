
#' Generate DNAString
#'
#' This function generates a random DNAString of a given length
#' by sampling from the given background model.
#'
#' @param len Integer length of the sequence
#' @inheritParams backgroundValid
#'
#' @return A DNAString object
#
#' @examples
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # generates a 1 kb sequence
#' motifcounter:::generateDNAString(1000,bg)
#'
#' @seealso \code{\link{generateDNAStringSet}}
generateDNAString=function(len,bg) {
    len=as.integer(len)
    backgroundValid(bg)
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
#' This function generates a DNAStringSet-object of given
#' individual sequence lengths
#' by sampling from the background model.
#'
#'
#' @inheritParams compoundPoissonDist
#' @inheritParams backgroundValid
#'
#' @return A DNAStringSet object
#
#' @examples
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load an order-1 background model
#' bg=readBackground(seqs,1)
#'
#' # generate sequences of various lengths
#' motifcounter:::generateDNAStringSet(10:50,bg)
#'
#' @seealso \code{\link{generateDNAStringSet}}
generateDNAStringSet=function(seqlen,bg) {
    backgroundValid(bg)
    seqs=c()
    for (i in 1:length(seqlen)) {
        seqs=c(seqs,generateDNAString(seqlen[i],bg))
    }

    return(Biostrings::DNAStringSet(seqs))
}

#' Empirical number of motif hits distribution
#'
#' This function repeatedly simulates random DNA sequences according to
#' the background model and
#' subsequently counts how many motif hits occur in them, 
#' which gives rise to the empirical distribution of the number of motif
#' hits.
#' This function is only used for benchmarking analysis.
#'
#' @inheritParams numMotifHits
#' @inheritParams compoundPoissonDist
#' @param nsim Integer number of random samples.
#' @return A List that contains
#' \describe{
#' \item{dist}{Empirical number of motif hits distribution}
#' }
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
#' simc=motifcounter:::simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=FALSE)
#'
#' simc=motifcounter:::simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=TRUE)
#'
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
simulateNumHitsDist=function(pfm,bg,seqlen, nsim, singlestranded=FALSE) {
    motifValid(pfm)
    backgroundValid(bg)
    if (length(nsim)<=0) {
        stop("nsim must be strictly positive")
    }

    freq=rep(0,10)
    for (i in 1:nsim) {
        seqs=generateDNAStringSet(seqlen,bg)
        nom=numMotifHits(seqs,pfm,bg,singlestranded)
        nom=sum(nom$numofhits)
        if (nom>=length(freq)) {
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

#' Empirical score distribution
#'
#' This function estimates the empirical score distribution
#' by simulating random DNA sequences based on the
#' background model. Subsequently, the random sequences are scanned
#' with the scoring measure which yields the score distribution.
#' This function is only used for benchmarking analysis.
#' 
#'
#' @inheritParams simulateNumHitsDist
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
#' # Load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # generate the simulated score distribution on
#' # sequences of length 1kb using 1000 samples
#' motifcounter:::scoreDistEmpirical(motif,bg,seqlen= 1000,nsim=1000)
#'
#' @seealso \code{\link{scoreDist}}
scoreDistEmpirical=function(pfm,bg,seqlen, nsim) {
  motifValid(pfm)
  backgroundValid(bg)
  seqs=generateDNAStringSet(rep(seqlen,nsim),bg)
  sh=scoreHistogram(seqs,pfm,bg)
  
  probs=sh$frequency/sum(sh$frequency)
  return(list(score=sh$score, probability=probs))
}
