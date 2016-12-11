#' Motif hit observations
#'
#' This function determines per-position motif hits in a given DNA sequence
#' @include score_wrapper.R
#' @inheritParams scoreSequence
#' @return List containing
#' \describe{
#' \item{fhits}{Per-position motif hits on the forward strand}
#' \item{rhits}{Per-position motif hits on the reverse strand}
#' }
#'
#' @examples
#'
#'
#' # Load DNA sequences
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seq=Biostrings::readDNAStringSet(seqfile)
#'
#' # Load a background model
#' bg=readBackground(seq,1)
#'
#' # Load a motif
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Determine the motif hits
#' motifHits(seq[[1]],motif,bg)
#'
#' @export
motifHits=function(seq,pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)

    sth=scoreThreshold(pfm,bg)
    scores=scoreSequence(seq,pfm,bg)
    fhits=rep(0,length(scores$fscores))
    rhits=rep(0,length(scores$rscores))
    fhits[scores$fscores>=sth$threshold]=1
    rhits[scores$rscores>=sth$threshold]=1

    return(list(fhits=fhits,rhits=rhits))
}

#' Motif hit profile across multiple sequences
#'
#' This function computes the per-position average motif hit 
#' profile across a set of fixed-length DNA sequences.
#'
#' @inheritParams scoreSequenceProfile
#'
#' @return List containing
#' \describe{
#' \item{fscores}{Per-position average forward strand motif hits}
#' \item{rscores}{Per-position average reverse strand motif hits}
#' }
#'
#' @examples
#'
#'
#'
#' # Load DNA sequences
#' seqfile=system.file("extdata","oct4_chipseq.fa", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' seqs=seqs[1:10]
#'
#' # Estimate a background model
#' bg=readBackground(seqs,1)
#'
#' # Load a motif
#' motiffile=system.file("extdata","x31.tab",package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the motif hit profile
#' motifHitProfile(seqs,motif,bg)
#'
#' @export
motifHitProfile=function(seqs,pfm,bg) {
    motifValid(pfm)
    backgroundValid(bg)
    if (class(seqs)!="DNAStringSet") {
        stop("seq must be a DNAString object")
    }
    if (any(lenSequences(seqs)!=length(seqs[[1]]))) {
        stop("all sequences must be equally long")
    }
    
    fhits=sapply(seqs, function(seq,pfm,bg) {
        motifHits(seq,pfm,bg)$fhits}, 
        pfm,bg)
    fhits=apply(fhits,1,mean)
    
    rhits=sapply(seqs, function(seq,pfm,bg) {
        motifHits(seq,pfm,bg)$rhits}, 
        pfm,bg)
    rhits=apply(rhits,1,mean)
    return (list(fhits=fhits,rhits=rhits))
}

#' Number of motif hits in a given DNA sequence
#'
#' This function counts the number of motif hits that
#' are found in the supplied DNA sequences.
#' It can be used to count motif hits on 
#' one or both strands, respectively.
#'
#'
#' @inheritParams scoreSequenceProfile
#' @param singlestranded Boolean that indicates whether a single strand or
#' both strands shall be scanned for motif hits. 
#' Default: singlestranded = FALSE.
#' @return A list containing
#' \describe{
#' \item{nseq}{Number of individual sequences}
#' \item{lseq}{Vector of individual sequence lengths}
#' \item{numofhits}{Vector of the number of hits in each individual sequence}
#' }
#' @examples
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' seq=Biostrings::readDNAStringSet(seqfile)
#'
#' # Set false positive probability
#' alpha=0.001
#' motifcounterOption(alpha)
#'
#' # Estimate an order-1 background model
#' bg=readBackground(seq,1)
#' # read PFM from file
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # scan the given sequence on both strands for the motif occurrences
#' noc=motifcounter:::numMotifHits(seq,motif,bg)
#' noc
#'
#' # scan the given sequence on a single strand for the motif occurences
#' noc=motifcounter:::numMotifHits(seq,motif,bg,singlestranded=TRUE)
#' noc
#'
numMotifHits=function(seqs, pfm, bg, singlestranded=FALSE) {
    motifValid(pfm)
    backgroundValid(bg)
    if (class(seqs)=="DNAStringSet") {
        # retrieve the number of motif hits
        x=lapply(seqs, function(seq,pfm,bg,singlestranded) {
            ret=motifHits(seq,pfm,bg)
            if (singlestranded==FALSE) {
                return(sum(ret[[1]]+ret[[2]]))
            } else {
                return(sum(ret[[1]]))
            }
        }, pfm,bg,singlestranded)

        noh=unlist(x)
        # retrieve the individual sequence lengths
        # sequences containing "N" or "n" are assigned length zero
        lseq=lenSequences(seqs)

        nseq=length(seqs)
    } else {
        stop("seqs must be a DNAStringSet object")
    }
    return (list(nseq=nseq, lseq=lseq,
        numofhits=noh))
}
