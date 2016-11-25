#' Estimates the background model
#' 
#' This function reads a DNA sequence from a given fasta file
#' and uses that sequence to estimate an order-m Markov model.
#' 
#' 
#' @param file Fasta-filename.
#' @param order Order of the Markov models that shall be used as the
#' background model. Default: order=1.
#'
#' @return None
#'
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
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
#' This function returns the parameters of the current background model
#' as an R list.
#' 
#' @return A list containing the stationary
#' and the transition probabilities, respectively:
#' \describe{
#' \item{stat}{Stationary distribution}
#' \item{trans}{Transition probabilities}
#' }
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' file=system.file("extdata","seq.fasta", package="mdist")
#' readBackground(file,1)
#' fetchBackground()
#' 
#' @export
fetchBackground=function() {
    stat=.Call("mdist_fetchStationBackground",PACKAGE="mdist");
    trans=.Call("mdist_fetchTransBackground",PACKAGE="mdist");
    return(list(stat=stat,trans=trans))
}



#' Prints the current background model
#' 
#' This function prints the currently loaded background model.
#' The function is primarily used for debugging and testing.
#' 
#' @return None
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
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
#' memory.
#' 
#' @return None
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
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
#' This function reads a DNA sequence from a given fasta file
#' and uses that sequence to estimate an order-m Markov model.
#' \strong{Note}: This function is only used for generating the random DNA
#' sequence that is used for computing the empirical
#' distribution. When using  \code{\link{compoundPoissonDist}},
#' \code{\link{combinatorialDist}} or \code{\link{motifEnrichmentTest}},
#' this function is not relevant. Instead, consult
#' \code{link{readBackground}}.
#' 
#' 
#' @param file Fasta-filename.
#' @param order Order of the Markov models that shall be used as the
#' background model. Default: order=1.
#' 
#' @return None
#' 
#' @seealso \code{link{readBackground}}
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
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
#' Similar to \code{link{printBackgroundForSampling}}, but
#' prints parameters that where acquired using 
#' \code{link{readBackgroundForSampling}}.
#' 
#' 
#' @return None
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackgroundForSampling(file=seqfile,1)
#' printBackgroundForSampling()
#' 
#' 
#' @seealso \code{link{printBackgroundForSampling}}
#' @export
printBackgroundForSampling=function() {
    dummy=.C("mdist_printBackgroundForSampling",PACKAGE="mdist");
}



#' Delete the current background model for sampling
#' 
#' Similar to \code{link{deleteBackgroundForSampling}}, but
#' deletes parameters that where acquired using 
#' \code{link{readBackgroundForSampling}}.
#' 
#' 
#' @return None
#' 
#' 
#' @examples
#' 
#' # Estimate first order Markov model based on the sequence provided
#' # in seq.fasta
#' 
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' readBackgroundForSampling(file=seqfile,1)
#' deleteBackgroundForSampling()
#' 
#' @seealso \code{link{deleteBackgroundForSampling}}
#' @export
deleteBackgroundForSampling=function() {
    dummy=.C("mdist_deleteBackgroundForSampling",PACKAGE="mdist")
}
