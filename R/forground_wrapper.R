#' Check valididity of PFM
#' 
#' This function checks if the PFM is valid. The function throws
#' an error if the R matrix does not represent a PFM.
#' 
#' @param pfm An R matrix that represents a position frequency matrix
#' @return None
#' 
#' @examples
#' 
#' motiffile=system.file("extdata","x1.tab", package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#' motifcounter:::motifValid(motif)
#' 
motifValid=function(pfm) {
    if (!is.matrix(pfm)) {
        stop("pfm must be a matrix")
    }
    #check if matrix has four rows
    if (nrow(pfm)!=4) {
        stop("pfm ncol must equal 4")
    }
    #check if all entries are positive
    if (!all(pfm>0)) {
        stop("pfm must be strictly positive.
            add a small pseudocount in case they are not")
    }
    #check if all columns sum to one
    if (!all(abs(1-apply(pfm,2,sum))<0.0000001)) {
        stop("all columns must sum to one, which seems not
            to be the case")
    }
}

#' Normalizes a PFM
#' 
#' This function normalizes a PFM and optionally
#' adds pseudo-evidence to each entry of the matrix.
#' 
#' @inheritParams motifValid
#' @param pseudo Small numeric pseudo-value that is added 
#' to each entry in the PFM in order to ensure strictly positive entries.
#' Default: pseudo = 0.01
#' @return None
#' 
#' @examples
#' 
#' 
#' motiffile=system.file("extdata","x1.tab", package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#' new_motif=normalizeMotif(motif)
#' 
#' @export
normalizeMotif=function(pfm,pseudo=0.01) {
    pfm=pfm+pseudo
    pfm=t(t(pfm)/apply(pfm,2,sum))
    return(pfm)
}

