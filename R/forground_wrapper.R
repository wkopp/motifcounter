#' Check valididity of PFM
#'
#' This function checks if the PFM is valid. The function throws
#' an error if the R matrix does not represent a PFM.
#'
#' @param pfm A position frequency matrix
#' @return None
#'
#' @examples
#'
#' motiffile=system.file("extdata","x1.tab", package="motifcounter")
#' motif=t(as.matrix(read.table(motiffile)))
#' motifValid(motif)
#'
#' @export
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

#' Add pseudo-count to PFM
#'
#' This function can be used to prevent the PFM from
#' containing zero-value entries.
#'
#' @param pfm A position frequency matrix
#' @param pseudo Small pseudo-value that is added to each entry in the PFM
#' @return A renormalized pfm
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
