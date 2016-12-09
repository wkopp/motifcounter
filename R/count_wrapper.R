#' Motif hit observations
#'
#' This function determines motif hits in a given DNA sequence
#'
#' @param seq A DNAString
#' @param pfm A position frequency matrix
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{fhits}{Vector of motif hits on the forward strand}
#' \item{rhits}{Vector of motif hits on the reverse strand}
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
#' bg=readBackground(seq,1)
#'
#' # Load the motif from the motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' scoreSequence(seq[[1]],motif,bg)
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
#' This function computes the average motif hit 
#' profile across a set of fixed-length DNA sequences.
#'
#' @param pfm A position frequency matrix
#' @param seqs DNAStringSet
#' @param bg A Background object
#' @return List containing
#' \describe{
#' \item{fscores}{Average forward strand motif hits}
#' \item{rscores}{Average reverse strand motif hits}
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
#' # Estimate the background model
#' bg=readBackground(seqs,1)
#'
#' # Load the motif from the motiffile
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
#' @param seqs A DNAString or a DNAStringSet object
#' @param pfm A position frequency matrix
#' @param bg A Background object
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits. Default: singlestranded=FALSE.
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
#' noc=numMotifHits(seq,motif,bg)
#' noc
#'
#' # scan the given sequence on a single strand for the motif occurences
#' noc=numMotifHits(seq,motif,bg,singlestranded=TRUE)
#' noc
#'
#' @export
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
        lseq=sapply(seqs, function(seq) {
            if (Biostrings::countPattern("N",seq)>0 ||
                Biostrings::countPattern("n",seq)>0) {
                return (0)
            } else {
                return(length(seq))
            }
        })

        nseq=length(seqs)
    } else if (class(seqs)=="DNAString") {
        nseq=1
        if (Biostrings::countPattern("N",seqs)>0 ||
            Biostrings::countPattern("n",seqs)>0) {
            lseq=0
            noh=0
        } else {
            lseq=length(seqs)
            ret=motifHits(seqs,pfm,bg)
            if (singlestranded==FALSE) {
                noh=sum(ret[[1]]+ret[[2]])
            } else {
                noh=sum(ret[[1]])
            }
        }
    } else {
        stop("seqs must be a DNAStringSet or DNAString object")
    }
    return (list(nseq=nseq, lseq=lseq,
        numofhits=noh))
}
