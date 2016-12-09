#' Enrichment of motif hits
#'
#' This function determines whether a given motif is enriched in a given
#' DNA sequences. Enrichment is tested by comparing the observed 
#' number of motif hits against the distribution of the number
#' of motif hits in random DNA sequences.
#' The function approximates the distribution of the number of motif
#' hits using either a 'compound Poisson approximation' 
#' or the 'combinatorial model'.
#'
#'
#' @param seqs DNAString or DNAStringSet
#' @param pfm A position frequency matrix
#' @param bg A Background object
#' @param singlestranded Boolean flag that indicates whether one or both
#'      strands are scanned for motif hits
#' @param method String that defines whether to use
#' the 'compound' Poisson approximation' or the 'combinatorial' model.
#' Default: method='compound'.
#'
#' @return Result list that contains
#' \describe{
#' \item{pvalue}{P-value for the enrichment test}
#' \item{fold}{Fold-enrichment with respect to the expected number of hits}
#' }
#' @examples
#'
#'
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' seqs=Biostrings::readDNAStringSet(seqfile)
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' alpha=0.001
#' gran=0.1
#' motifcounterOption(alpha, gran)
#'
#' # estimate the background model
#' bg=readBackground(seqs,1)
#'
#' # load a motif
#' motif=t(as.matrix(read.table(motiffile)))
#'
#' ### 1 ) Compute the distribution for scanning a *single* DNA strand
#' # based on the 'Compound Poisson model'
#'
#' result=motifEnrichment(seqs,motif,bg,
#'             singlestranded=TRUE,method="compound")
#'
#' ### 2 ) Compute the distribution for scanning *both* DNA strand
#' # based on the 'Compound Poisson model'
#'
#' result=motifEnrichment(seqs,motif, bg, method="compound")
#'
#' ### 3 ) Compute the distribution for scanning *both* DNA strand
#' # based on the *combinatorial model*
#'
#' result=motifEnrichment(seqs,motif, bg,singlestranded=FALSE,
#'             method="combinatorial")
#'
#' @seealso \code{\link{compoundPoissonDist}}, \code{\link{combinatorialDist}}
#' @export
motifEnrichment=function(seqs, pfm,bg,
    singlestranded=FALSE,method="compound") {
    motifValid(pfm)
    backgroundValid(bg)
    #compute overlapping hit probs
    overlap=probOverlapHit(pfm,bg,singlestranded)

    # detemine the number of motif hits
    observations=numMotifHits(seqs,pfm,bg,singlestranded)
    if (method=="compound") {
        dist=compoundPoissonDist(observations$lseq, overlap)
    } else if (method=="combinatorial") {
        dist=combinatorialDist(observations$lseq, overlap)
    } else {
        stop("method must be 'compound' or 'combinatorial'")
    }
    p=sum(dist$dist[(sum(observations$numofhits)+1):length(dist$dist)])

    return (list(pvalue=p,fold=sum(observations$numofhits)/
            sum(dist$dist*seq(0,length(dist$dist)-1))))
}
