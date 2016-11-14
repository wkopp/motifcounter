#' Read and load the PFM
#' 
#' The function reads a position frequency matrix (PFM) for scanning a DNA
#' sequence. The PFM can be loaded directly in transfac format, in tab format
#' or from an R matrix. To load a PFM from a file directly, the files must end
#' in '.transfac' or '.tab' to load from transfac or tab format, respectively.
#' 
#' 
#' @param pwm Character string of the name of the file which contains the PFM
#' or an R matrix.
#' @param pseudocount Numeric value which is added to each element of the
#' loaded PFM in order to avoid zero probability elements. The columns of the
#' matrix are renormalized automatically.
#' @examples
#' 
#' 
#' library(mdist)
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' # load a motif in tab format
#' readMotif(pwmfile, 0.01)
#' 
#' # fetch the motif to an R matrix
#' pwm=motif2matrix()
#' 
#' #reload the motif again from the R matrix
#' readMotif(pwm, 0.01)
#' 
#' 
#' @export
readMotif=function(pwm, pseudocount=0.01) {
  if (is.matrix(pwm)) {
    if (any(pwm<=0)) {
        stop("All entries in the matrix must be greater than zero")
    }
    if (nrow(pwm)!=4) {
        stop("The number of rows must be 4, representing the number nucleotides.")
    }
    pwm=pwm+pseudocount
    pwm=pwm/apply(pwm,2,sum)
    dummy=.C("mdist_loadmotif", as.numeric(pwm), nrow(pwm), ncol(pwm),PACKAGE="mdist")
  } else if (is.character(pwm)) {
    sdummy=.C("mdist_motiffromfile", as.character(pwm),
        as.numeric(pseudocount),PACKAGE="mdist")
  } else {
      stop("pwm must be a filename pointing that contains a PFM or a PFM matrix")
  }
}



#' Deletes the current motif
#' 
#' This function unloads the current motif and frees the allocated memory
#' 
#' 
#' @examples
#' 
#' 
#' library(mdist)
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' readMotif(pwmfile,0.01)
#' deleteMotif()
#' 
#' @export
deleteMotif=function() {
  dummy=.C("mdist_deleteMotif",PACKAGE="mdist")
}



#' Get the motif in an R matrix
#' 
#' This function returns the current PFM as an R matrix.
#' 
#' 
#' @examples
#' 
#' 
#' library(mdist)
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' readMotif(pwmfile)
#' x=motif2matrix()
#' 
#' @export
motif2matrix=function() {
	m=.Call("mdist_fetchMotif",PACKAGE="mdist");
	return (m)
}

#' Length of the currently loaded PFM
#' 
#' This function returns the length of the currently loaded PFM.
#' 
#' 
#' @examples
#' 
#' 
#' library(mdist)
#' pwmfile=system.file("extdata","x31.tab", package="mdist")
#' readMotif(pwmfile)
#' motifLength()
#' 
#' @export
motifLength=function() {
  x=integer(1);
  res=.C("mdist_motiflength",x,PACKAGE="mdist")
  return (res[[1]])
}

