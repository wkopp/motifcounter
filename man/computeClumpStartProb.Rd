% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markovmodel.R
\name{computeClumpStartProb}
\alias{computeClumpStartProb}
\title{Computes the Clump start probability based on a Markov model}
\usage{
computeClumpStartProb(overlap)
}
\arguments{
\item{overlap}{An Overlap object.}
}
\value{
Clump start probability 'tau'
}
\description{
This function leverages a Markov model in order to
determine the clump start probability.
The computation depends on the selected false positive probability 
for calling motif matches 'alpha' and the pre-determined
overlapping match probabilities 'beta'.
}
\details{
The general idea of the method relies on the fact that
for the stationary distribution of the Markov model,
motif matches must be observed with probability 'alpha'.
Hence, the clump start probability 'tau' is optimized
to achieve that goal.

The R interface is only used for the purpose of testing
the correctness of the model.
}
\examples{

# Load sequences
seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
seqs = Biostrings::readDNAStringSet(seqfile)

# Load motif
motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
motif = t(as.matrix(read.table(motiffile)))

# Load background model
bg = readBackground(seqs, 1)

# Compute overlap probabilities
op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)


# Computes the clump start probability
dist = motifcounter:::computeClumpStartProb(op)

}
\seealso{
\code{\link{compoundPoissonDist}}

\code{\link{numMotifHits}}

\code{\link{probOverlapHit}}
}
