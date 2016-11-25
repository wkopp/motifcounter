#' Read and load a PFM
#' 
#' The function reads a position frequency matrix (PFM) for scanning a DNA
#' sequence. The PFM can be loaded directly in transfac format, in tab format
#' or from an R matrix. To load a PFM from a file directly, the files must end
#' in '.transfac' or '.tab' to load from transfac or tab format, respectively,
#' and each file must contain only one PFM.
#' 
#' 
#' @param motif Filename to load the PFM from
#' or PFM provided as an R matrix.
#' @param pseudocount Numeric value which is added to each element of the
#' PFM in order to avoid zero probability elements. The columns of the
#' matrix are renormalized automatically.
#' 
#' @return None
#' @examples
#' 
#' 
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' # load a motif in tab format
#' readMotif(motiffile, 0.01)
#' 
#' # fetch the motif to an R matrix
#' motif=motif2matrix()
#' 
#' #reload the motif again from the R matrix
#' readMotif(motif, 0.01)
#' 
#' 
#' @export
readMotif=function(motif, pseudocount=0.01) {
    if (is.matrix(motif)) {
        if (any(motif<=0)) {
            stop("All entries in the matrix must be greater than zero")
        }
        if (nrow(motif)!=4) {
            stop("The number of rows must be 4, 
                representing the number nucleotides.")
        }
        motif=motif+pseudocount
        motif=motif/apply(motif,2,sum)
        dummy=.C("motifcounter_loadmotif", as.numeric(motif), 
                nrow(motif), ncol(motif),PACKAGE="motifcounter")
    } else if (is.character(motif)) {
        sdummy=.C("motifcounter_motiffromfile", as.character(motif),
                as.numeric(pseudocount),PACKAGE="motifcounter")
    } else {
        stop("motif must be a filename pointing 
                that contains a PFM or a PFM matrix")
    }
}



#' Deletes the current motif
#' 
#' This function unloads the current motif and frees the allocated memory
#' 
#' @return None
#' 
#' @examples
#' 
#' 
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' readMotif(motiffile,0.01)
#' deleteMotif()
#' 
#' @export
deleteMotif=function() {
    dummy=.C("motifcounter_deleteMotif",PACKAGE="motifcounter")
}



#' Returns the current PFM as an R matrix
#' 
#' This function returns the current PFM as an R matrix.
#' 
#' @return  A position frequency matrix
#' 
#' @examples
#' 
#' 
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' readMotif(motiffile)
#' x=motif2matrix()
#' 
#' @export
motif2matrix=function() {
    m=.Call("motifcounter_fetchMotif",PACKAGE="motifcounter");
    return (m)
}

#' Length of the currently loaded PFM
#' 
#' This function returns the length of the currently loaded PFM.
#' 
#' @return Length of the motif
#' 
#' @examples
#' 
#' 
#' motiffile=system.file("extdata","x31.tab", package="motifcounter")
#' readMotif(motiffile)
#' motifLength()
#' 
#' @export
motifLength=function() {
    x=integer(1);
    res=.C("motifcounter_motiflength",x,PACKAGE="motifcounter")
    return (res[[1]])
}

