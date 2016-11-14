#' Estimates the background model
#' 
#' This function reads a DNA sequence in fasta format and estimates and order m
#' Markov model.
#' 
#' 
#' @param file Filename where the sequence is stored in. The file must be in
#' fasta format.
#' @param order Defines the order of the Markov models used for the background.
#' Default: order=1.
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' file=system.file("extdata","seq.fasta", package="mdist")
#' readBackground(file,1)
#' 
#' @export
readBackground=function(file, order=1) {
   nseq=numSequences(file)
   lseq=lenSequences(file)
   dummy=.C("mdist_makebg", as.character(file), as.integer(order),
    				 as.integer(nseq),as.integer(lseq),PACKAGE="mdist")
}

#' Fetch background model into R vectors
#' 
#' This function fetches the current stationary distribution
#' and transition probabilities into R two vectors that are
#' contained in list. This function is only used for debugging 
#' and testing purposes.
#' 
#' @return A list containing a vector containing the stationary
#'              and the transition probabilities, respectively.
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' file=system.file("extdata","seq.fasta", package="mdist")
#' readBackground(file,1)
#' fetchBackground()
#' 
#' @export
fetchBackground=function() {
    stat=.Call("mdist_fetchStationBackground",PACKAGE="mdist");
    trans=.Call("mdist_fetchTransBackground",PACKAGE="mdist");
    return(list(stat,trans))
}



#' Prints the current background model
#' 
#' This function prints the currently loaded background model.
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackground(file=seqfile,1)
#' printBackground()
#' 
#' @export
printBackground=function() {
  dummy=.C("mdist_printBackground",PACKAGE="mdist");
}



#' Delete background model
#' 
#' This function unloads the current background model and frees its allocated
#' memeory.
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackground(file=seqfile,1)
#' deleteBackground()
#' 
#' 
#' @export
deleteBackground=function() {
  dummy=.C("mdist_deleteBackground",PACKAGE="mdist")
}



#' Estimates the background model for sampling
#' 
#' This function reads a DNA sequence in fasta format and estimates and order m
#' Markov model. This function is only used when simulated DNA sequences are
#' scanned by the model.  The compound Poisson model and the combinatorial
#' model do not use this backgroud model.
#' 
#' 
#' @param file Name of the sequence file in fasta format.
#' @param order Defines the order of the Markov models used for the background.
#' Default: order=1.
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' file=system.file("extdata","seq.fasta", package="mdist")
#' readBackgroundForSampling(file,1)
#' 
#' 
#' @export
readBackgroundForSampling=function(file, order=1) {
	 nseq=numSequences(file)
   lseq=lenSequences(file)
  dummy=.C("mdist_makebgForSampling", as.character(file), as.integer(order),
  				 as.integer(nseq),as.integer(lseq),PACKAGE="mdist")
}



#' Prints the current background model for sampling
#' 
#' This function prints the currently loaded background model that is employed
#' during scanning the sampled DNA sequence. It is not used for the compound
#' Poisson model and the combinatorial model
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackgroundForSampling(file=seqfile,1)
#' printBackgroundForSampling()
#' 
#' 
#' @export
printBackgroundForSampling=function() {
  dummy=.C("mdist_printBackgroundForSampling",PACKAGE="mdist");
}



#' Delete the current background model for sampling
#' 
#' This function unloads the current background model for sampling.
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' library(mdist)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackgroundForSampling(file=seqfile,1)
#' deleteBackgroundForSampling()
#' 
#' 
#' @export
deleteBackgroundForSampling=function() {
  dummy=.C("mdist_deleteBackgroundForSampling",PACKAGE="mdist")
}
