% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_wrapper.R
\name{scoreThreshold}
\alias{scoreThreshold}
\title{Score threshold}
\usage{
scoreThreshold(pfm, bg)
}
\arguments{
\item{pfm}{An R matrix that represents a position frequency matrix}

\item{bg}{A Background object}
}
\value{
List containing
\describe{
\item{threshold}{Score threshold}
\item{alpha}{False positive probability}
}
}
\description{
This function computes the score threshold for a desired
false positive probability `alpha`.
}
\details{
Note that the returned alpha usually differs slightly
from the one that is prescribed using
\code{\link{motifcounterOptions}}, because
of the discrete nature of the sequences.
}
\examples{

# Load sequences
seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
seqs = Biostrings::readDNAStringSet(seqfile)

# Load background
bg = readBackground(seqs, 1)

# Load motif
motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
motif = t(as.matrix(read.table(motiffile)))

# Compute the score threshold
motifcounter:::scoreThreshold(motif, bg)

}
