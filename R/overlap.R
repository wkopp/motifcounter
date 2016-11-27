#' Overlapping motif hit probabilities
#' 
#' This function computes the self-overlapping probabilites of the motif for
#' the current PFM and background based on two-dimensional score distributions.
#' Subsequently, the overlapping hit probabilities are corrected for motif hits
#' that would arise in between the shifted positions.
#' 
#' @param pfm A position frequency matrix
#' @param singlestranded Boolian which defines whether the overlapping hit
#' probabilities shall be computed with respect to scanning both DNA strands or
#' only one strand.  Default: Both strands are scanned for DNA sequences.
#' @return A list containing various overlapping hit probabilities.
#' The list contains the following entries
#'     \describe{
#'     \item{alpha}{False positive motif hit probability}
#'     \item{beta}{Vector of overlapping hit probability for hits 
#'             occurring on the same strand. Each element corresponds to
#'             relative distance of the starting positions between the
#'             motif hits.}
#'     \item{beta3p}{Vector of overlapping hit probability for a 
#'             forward strand hit that is followed by a reverse strand
#'             hit. Each element corresponds to
#'             relative distance of the starting positions between the
#'             motif hits.}
#'     \item{beta5p}{Vector of overlapping hit probability for a
#'             reverse strand hit that is followed by a forward strand
#'             hit. Each element corresponds to
#'             relative distance of the starting positions between the
#'             motif hits.}
#'     \item{gamma}{Vector of overlapping hit probabilities across
#'             all configurations. In contrast to beta, beta3p and beta5p,
#'             gamma is not corrected for intermediate motif hit events.
#'             only used for the compound Poisson variant according to Pape
#'             et al. 2008.}
#'     \item{singlestranded}{returns to singlestranded flag}
#'     }
#'              
#' @examples
#' 
#' 
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' alpha=0.001
#' gran=0.1
#' motifcounterOption(alpha, gran)
#' 
#' # estimate background model from seqfile
#' readBackground(seqfile,1)
#' 
#' # load motif model from motiffile
#' motif=t(as.matrix(read.table(motiffile)))
#' 
#' # compute the overlap probabilities for scanning both DNA strands
#' op=probOverlapHit(motif,singlestranded=FALSE)
#' 
#' # compute the overlap probabilities for scanning a single DNA strand
#' op=probOverlapHit(motif,singlestranded=TRUE)
#' 
#' @export
probOverlapHit=function(pfm,singlestranded=FALSE) {
    #check if pfm is a matrix
    motifValid(pfm)
    alpha=numeric(1)
    beta=numeric(ncol(pfm))
    beta3p=numeric(ncol(pfm))
    beta5p=numeric(ncol(pfm))
    gamma=numeric(3*ncol(pfm))
    if (singlestranded==TRUE) {
        res=.C("motifcounter_overlapSingleStranded", 
            as.numeric(pfm),nrow(pfm), ncol(pfm), 
            alpha, beta, beta3p, beta5p, gamma,PACKAGE="motifcounter")
    } else {
        res=.C("motifcounter_overlap", 
            as.numeric(pfm),nrow(pfm), ncol(pfm), 
            alpha, beta, beta3p, beta5p, gamma,PACKAGE="motifcounter")
    }
    return (list(alpha=res[[4]],beta=res[[5]],beta3p=res[[6]],beta5p=res[[7]],
        gamma=res[[8]], singlestranded=singlestranded))
}

