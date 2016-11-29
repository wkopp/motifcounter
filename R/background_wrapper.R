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
#' file=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackground(file,1)
#'
#' @export
readBackground=function(file, order=1) {
    seq=Biostrings::readDNAStringSet(file)
    nseq=numSequences(seq)
    lseq=lenSequences(seq)
    dummy=.C("motifcounter_makebg", as.character(file), as.integer(order),
        as.integer(nseq),as.integer(lseq),PACKAGE="motifcounter")
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
#' file=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackground(file,1)
#' fetchBackground()
#'
#' @export
fetchBackground=function() {
    stat=.Call("motifcounter_fetchStationBackground",PACKAGE="motifcounter");
    trans=.Call("motifcounter_fetchTransBackground",PACKAGE="motifcounter");
    return(list(stat=stat,trans=trans))
}

#' Get background order
#'
#' This function returns the background model order
#'
#' @return Order of the background model
#'
#'
bgOrderForSampling=function() {
    order=integer(1)
    order=.C("motifcounter_bgorder",as.integer(order),
        PACKAGE="motifcounter");
    return(order[[1]])
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
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackground(file=seqfile,1)
#' printBackground()
#'
#' @export
printBackground=function() {
    dummy=.C("motifcounter_printBackground",PACKAGE="motifcounter");
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
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackground(file=seqfile,1)
#' deleteBackground()
#'
#'
#' @export
deleteBackground=function() {
    dummy=.C("motifcounter_deleteBackground",PACKAGE="motifcounter")
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
#' file=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackgroundForSampling(file,1)
#'
#'
#' @export
readBackgroundForSampling=function(file, order=1) {
    seq=Biostrings::readDNAStringSet(file)
    nseq=numSequences(seq)
    lseq=lenSequences(seq)
    dummy=.C("motifcounter_makebgForSampling",
        as.character(file), as.integer(order),
        as.integer(nseq),as.integer(lseq),PACKAGE="motifcounter")
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
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackgroundForSampling(file=seqfile,1)
#' printBackgroundForSampling()
#'
#'
#' @seealso \code{link{printBackgroundForSampling}}
#' @export
printBackgroundForSampling=function() {
    dummy=.C("motifcounter_printBackgroundForSampling",PACKAGE="motifcounter");
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
#' seqfile=system.file("extdata","seq.fasta", package="motifcounter")
#' readBackgroundForSampling(file=seqfile,1)
#' deleteBackgroundForSampling()
#'
#' @seealso \code{link{deleteBackgroundForSampling}}
#' @export
deleteBackgroundForSampling=function() {
    dummy=.C("motifcounter_deleteBackgroundForSampling",PACKAGE="motifcounter")
}
