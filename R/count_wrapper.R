#' Motif hit observations
#'
#' This function computes the observed motif hits in a given DNA sequence
#'
#' @param pfm A position frequency matrix
#' @param seq DNAString
#' @return List containing
#' \describe{
#' \item{fhits}{Motif hits on the forward strand}
#' \item{rhits}{Motif hits on the reverse strand}
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
motifHits=function(pfm,seq) {
    motifValid(pfm)
    sth=scoreThreshold(pfm)
    scores=scoreSequence(pfm,seq)
    fhits=rep(0,length(scores$fscores))
    rhits=rep(0,length(scores$rscores))
    fhits[scores$fscores>=sth$threshold]=1
    rhits[scores$rscores>=sth$threshold]=1

    return(list(fhits=fhits,rhits=rhits))
}


#' Number of motif hits in a given DNA sequence
#'
#' This function scans the DNA sequences contained in the fasta file
#' and counts the number of motif hits using the score threshold
#' that is associated with the false positive probability 'alpha'
#' (see \code{\link{motifcounterOption}}. The function can be used
#' to count motif hits on one or both strands, respectively.
#'
#'
#' @param pfm A position frequency matrix
#' @param seqs A DNAString or a DNAStringSet object
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits
#' @return A list containing
#' \describe{
#' \item{nseq}{Number of individual sequences}
#' \item{lseq}{Vector of individual sequence lengths}
#' \item{numofhits}{Vector of the number of hits in each individual sequence}
#' }
#' @examples
#'
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' seq=Biostrings::readDNAStringSet(seqfile)
#'
#' # Set false positive probability
#' alpha=0.001
#' motifcounterOption(alpha)
#'
#' # Estimate order-1 background model
#' readBackground(seqfile,1)
#' # read PFM from file
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' # scan the given sequence on both strands for the motif occurances
#' noc=numMotifHits(motif,seq)
#' noc
#'
#' # scan the given sequence on a single strand for the motif occurances
#' noc=numMotifHits(motif,seq,singlestranded=TRUE)
#' noc
#'
#' @export
numMotifHits=function(pfm, seqs, singlestranded=FALSE) {
    motifValid(pfm)
    if (class(seqs)=="DNAStringSet") {
        # retrieve the number of motif hits
        x=lapply(seqs, function(seq,singlestranded) {
            ret=motifHits(pfm,seq)
            if (singlestranded==FALSE) {
                return(sum(ret[[1]]+ret[[2]]))
            } else {
                return(sum(ret[[1]]))
            }
        }, singlestranded)

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
        if (Biostrings::countPattern("N",seq)>0 ||
            Biostrings::countPattern("n",seq)>0) {
            lseq=0
            noh=0
        } else {
            lseq=length(seqs)
            ret=motifHits(pfm,seqs)
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
