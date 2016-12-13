#' Generate DNAString
#'
#' This function generates a random DNAString of a given length
#' by sampling from the background model.
#'
#' @param len Integer length of the sequence
#' @inheritParams backgroundValid
#'
#' @return A DNAString object
#
#' @examples
#'
#' # Load sequences
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg=readBackground(seqs,1)
#'
#' # Generate a 1 kb random sequence
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
#' This function generates a DNAStringSet-object of the
#' given individual sequence lengths
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
#' # Load sequences
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg=readBackground(seqs,1)
#'
#' # Generate random sequences of various lengths
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
#' subsequently counts how many motif hits occur in them.
#' Thus, this function gives rise to the empirical 
#' distribution of the number of motif hits.
#' This function is only used for benchmarking analysis.
#'
#' @inheritParams numMotifHits
#' @inheritParams compoundPoissonDist
#' @param nsim Integer number of random samples.
#' @return A List that contains
#' \describe{
#' \item{dist}{Empirical distribution of the number of motif hits}
#' }
#' @examples
#'
#'
#' # Load sequences
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg=readBackground(seqs,1)
#'
#' # Load motif
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' # Study the counts in one sequence of length 150 bp
#' seqlen=rep(150,1)
#' 
#' # Compute empirical distribution of the number of motif hits
#' # by scanning both strands using 100 samples
#' simc=motifcounter:::simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=FALSE)
#'
#' # Compute empirical distribution of the number of motif hits
#' # by scanning a single strand using 100 samples
#' simc=motifcounter:::simulateNumHitsDist(motif, bg, seqlen, nsim=100,singlestranded=TRUE)
#'
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
simulateNumHitsDist=function(pfm,bg,seqlen, nsim, singlestranded=FALSE) {
    motifValid(pfm)
    backgroundValid(bg)
    motifAndBackgroundValid(pfm,bg)
    stopifnot(nsim>0)

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
#' on a set of randomly generated DNA sequences based on the
#' background model.
#' This function is only used for benchmarking analysis.
#' 
#'
#' @inheritParams simulateNumHitsDist
#' @return List containing
#' \describe{
#' \item{scores}{Vector of scores}
#' \item{dist}{Score distribution}
#' }
#' @examples
#'
#' # Load sequences
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg=readBackground(seqs,1)
#'
#' # Load motif
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#'
#' # Compoute the empirical score distribution in
#' # sequences of length 1kb using 1000 samples
#' motifcounter:::scoreDistEmpirical(motif,bg,seqlen= 1000,nsim=1000)
#'
#' @seealso \code{\link{scoreDist}}
scoreDistEmpirical=function(pfm,bg,seqlen, nsim) {
  motifValid(pfm)
  backgroundValid(bg)
  motifAndBackgroundValid(pfm,bg)
  
  seqs=generateDNAStringSet(rep(seqlen,nsim),bg)
  sh=scoreHistogram(seqs,pfm,bg)
  
  probs=sh$frequency/sum(sh$frequency)
  return(list(scores=sh$score, dist=probs))
}
