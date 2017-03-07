#' Score distribution
#'
#' This function computes the score distribution for the given PFM and
#' background. The Score distribution is computed based on an efficient
#' dynamic programming algorithm.
#'
#' @inheritParams motifValid
#' @template templates
#' @return List that contains
#' \describe{
#' \item{scores}{Vector of scores}
#' \item{dist}{Score distribution}
#' }
#'
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' dp = scoreDist(motif, bg)
#'
#' @export
scoreDist = function(pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    scores = .Call(
        motifcounter_scorerange,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )

    dist = .Call(
        motifcounter_scoredist,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )
    return(list(scores = scores, dist = dist))
}


#' Score distribution
#'
#' This function computes the score distribution for a given PFM and
#' a background model.
#'
#' The result of this function is identical to \code{\link{scoreDist}},
#' however, the method employs a less efficient algorithm that
#' enumerates all DNA sequences of the length of the motif.
#' This function is only used for debugging and testing purposes
#' and might require substantial computational
#' resources for long motifs.
#'
#' @inheritParams scoreDist
#' @return List containing
#' \describe{
#' \item{scores}{Vector of scores}
#' \item{dist}{Score distribution}
#' }
#'
#' @seealso \code{\link{scoreDist}}
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score distribution
#' dp = motifcounter:::scoreDistBf(motif, bg)
#'
scoreDistBf = function(pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    scores = .Call(
        motifcounter_scorerange,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )

    dist = .Call(
        motifcounter_scoredist_bf,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )
    return(list(scores = scores, dist = dist))
}

#' Score strand
#'
#' This function computes the per-position
#' score in a given DNA strand.
#'
#' The function returns the per-position scores
#' for the given strand. If the sequence is too short,
#' it contains an empty vector.
#'
#' @inheritParams scoreDist
#' @param seq A DNAString object
#' @return
#' \describe{
#' \item{scores}{Vector of scores on the given strand}
#' }
#'
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the per-position and per-strand scores
#' motifcounter:::scoreStrand(seqs[[1]], motif, bg)
#'
scoreStrand = function(seq, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    # Check class
    stopifnot(is(seq, "DNAString"))

    scores = .Call(
        motifcounter_scoresequence,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        toString(seq),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )
    return(as.numeric(scores))
}

#' Score observations
#'
#' This function computes the per-position and per-strand
#' score in a given DNA sequence.
#'
#' @inheritParams scoreDist
#' @param seq A DNAString object
#' @return List containing
#' \describe{
#' \item{fscores}{Vector of scores on the forward strand}
#' \item{rscores}{Vector of scores on the reverse strand}
#' }
#'
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the per-position and per-strand scores
#' scoreSequence(seqs[[1]], motif, bg)
#'
#' @export
scoreSequence = function(seq, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    # Check class
    stopifnot(is(seq, "DNAString"))

    fscores = scoreStrand(seq, pfm, bg)
    rscores = scoreStrand(seq, revcompMotif(pfm), bg)
    return(list(fscores = fscores, rscores = rscores))
}

#' Score profile across multiple sequences
#'
#' This function computes the per-position and per-strand
#' average score profiles across a set of DNA sequences.
#' It can be used to reveal positional constraints
#' of TFBSs.
#'
#' @inheritParams scoreDist
#' @param seqs A DNAStringSet object
#'
#' @return List containing
#' \describe{
#' \item{fscores}{Vector of per-position average forward strand scores}
#' \item{rscores}{Vector of per-position average reverse strand scores}
#' }
#'
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score profile
#' scoreProfile(seqs, motif, bg)
#'
#' @export
scoreProfile = function(seqs, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    stopifnot (is(seqs, "DNAStringSet"))
    
    if (any(lenSequences(seqs) != lenSequences(seqs)[1])) {
        stop(paste(strwrap("All DNAStrings in 'seqs' must be of equal length.
            Please trim the sequences accordingly."), collapse = "\n"))
    }
    slen = lenSequences(seqs[1])

    if (slen <= ncol(pfm) - 1) {
        return (list(fscores = integer(0), rscores = integer(0)))
    }

    fscores = vapply(seqs, function(seq, pfm, bg) {
        s = scoreStrand(seq, pfm, bg)
    },
    FUN.VALUE = numeric(slen - ncol(pfm) + 1),
    pfm, bg)
    fscores = rowMeans(as.matrix(fscores))

    rscores = vapply(seqs, function(seq, pfm, bg) {
        s = scoreStrand(seq, revcompMotif(pfm), bg)
    },
    FUN.VALUE = numeric(slen - ncol(pfm) + 1),
    pfm, bg)
    rscores = rowMeans(as.matrix(rscores))
    return (list(
        fscores = as.vector(fscores),
        rscores = as.vector(rscores)
    ))
}

#' Score histogram on a single sequence
#'
#' This function computes the empirical score
#' distribution by normalizing the observed score histogram
#' for a given sequence.
#'
#'
#' @inheritParams scoreSequence
#' @return List containing
#' \describe{
#' \item{scores}{Vector of scores}
#' \item{dist}{Score distribution}
#' }
#'
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the per-position and per-strand scores
#' motifcounter:::scoreHistogramSingleSeq(seqs[[1]], motif, bg)
#'
scoreHistogramSingleSeq = function(seq, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    stopifnot(is(seq, "DNAString"))
    motifAndBackgroundValid(pfm, bg)

    scores = .Call(
        motifcounter_scorerange,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )

    dist = .Call(
        motifcounter_scorehistogram,
        as.numeric(pfm),
        nrow(pfm),
        ncol(pfm),
        toString(seq),
        bg@station,
        bg@trans,
        as.integer(bg@order)
    )
    result = list(scores = scores, dist = dist)
    
    return(result)
}


#' Score histogram
#'
#' This function computes the empirical score
#' distribution for a given set of DNA sequences.
#'
#' It can be used to compare the empirical score
#' distribution against the theoretical one (see \code{\link{scoreDist}}).
#'
#' @inheritParams scoreProfile
#' @return List containing
#' \describe{
#' \item{scores}{Vector of scores}
#' \item{dist}{Score distribution}
#' }
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the empirical score histogram
#' scoreHistogram(seqs, motif, bg)
#'
#' @seealso \code{\link{scoreDist}}
#' @export
scoreHistogram = function(seqs, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    stopifnot(is(seqs, "DNAStringSet"))

    # First, extract the score range
    his = lapply(seqs[1], scoreHistogramSingleSeq, pfm, bg)
    nseq = length(his)
    scores = his[[1]]$scores

    # Then extract the histogram
    his = vapply(seqs, function(seq, pfm, bg) {
        scoreHistogramSingleSeq(seq, pfm, bg)$dist
    }, FUN.VALUE = numeric(length(scores)), pfm, bg)
    
    freq = rowSums(his)
    result = list(scores = scores, dist = freq)
    
    return(result)
}

#' Score threshold
#'
#' This function computes the score threshold for a desired
#' false positive probability `alpha`.
#'
#' Note that the returned alpha usually differs slightly
#' from the one that is prescribed using
#' \code{\link{motifcounterOptions}}, because
#' of the discrete nature of the sequences.
#'
#' @inheritParams scoreDist
#' @return List containing
#' \describe{
#' \item{threshold}{Score threshold}
#' \item{alpha}{False positive probability}
#' }
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the score threshold
#' motifcounter:::scoreThreshold(motif, bg)
#'
scoreThreshold = function(pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    
    scoredist = scoreDist(pfm, bg)
    
    # find quantile
    ind = which(1 - cumsum(scoredist$dist) <= sigLevel())
    if (length(ind) <= 1) {
        stop(paste(strwrap(
            "The significance level 'alpha' is too stringent for the given 'pfm',
            which renders motif matches impossible to occur.
            Use 'motifcounterOptions' to prescribe a less stringent
            value for 'alpha'."), collapse = "\n"))
    }
    ind = tail(ind, -1)
    alpha = sum(scoredist$dist[ind])
    
    ind = min(ind)
    threshold = scoredist$scores[ind]
    
    return(list(threshold = threshold, alpha = alpha))
}


