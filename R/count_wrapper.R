#' Motif hit observations
#'
#' This function determines per-position motif hits in a given DNA sequence.
#' @include score_wrapper.R
#' @inheritParams scoreSequence
#' @return List containing
#' \describe{
#' \item{fhits}{Per-position motif hits on the forward strand}
#' \item{rhits}{Per-position motif hits on the reverse strand}
#' }
#'
#' @examples
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seq = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seq, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Determine the motif hits
#' motifHits(seq[[1]], motif, bg)
#'
#' @export
motifHits = function(seq, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    sth = scoreThreshold(pfm, bg)
    scores = scoreSequence(seq, pfm, bg)
    fhits = integer(length(scores$fscores))
    rhits = integer(length(scores$rscores))
    fhits[scores$fscores >= sth$threshold] = 1
    rhits[scores$rscores >= sth$threshold] = 1

    return(list(fhits = as.integer(fhits), rhits = as.integer(rhits)))
}

#' Motif hit profile across multiple sequences
#'
#' This function computes the per-position average motif hit
#' profile across a set of fixed-length DNA sequences.
#' It can be used to reveal positional constraints
#' of TFBSs.
#'
#' @inheritParams scoreProfile
#'
#' @return List containing
#' \describe{
#' \item{fscores}{Per-position average forward strand motif hits}
#' \item{rscores}{Per-position average reverse strand motif hits}
#' }
#'
#' @examples
#'
#'
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#' seqs = seqs[1:10]
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Load motif
#' motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
#' motif = t(as.matrix(read.table(motiffile)))
#'
#' # Compute the motif hit profile
#' motifHitProfile(seqs, motif, bg)
#'
#' @export
motifHitProfile = function(seqs, pfm, bg) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    stopifnot(class(seqs) == "DNAStringSet")

    if (any(lenSequences(seqs) != lenSequences(seqs)[1])) {
        stop(paste(strwrap("All DNAStrings in 'seqs' must be of equal length.
            Please trim the sequences such that they are equally long."), 
            collapse = "\n"))
    }
    slen = lenSequences(seqs[1])
    if (slen <= ncol(pfm) - 1) {
        return (list(fhits = integer(0), rhits = integer(0)))
    }

    fhits = vapply(seqs, function(seq, pfm, bg) {
        return(motifHits(seq, pfm, bg)$fhits)
    }, FUN.VALUE = integer(slen - ncol(pfm) + 1), pfm, bg)

    fhits = rowMeans(as.matrix(fhits))

    rhits = vapply(seqs, function(seq, pfm, bg) {
        mh = motifHits(seq, pfm, bg)$rhits
    }, FUN.VALUE = integer(slen - ncol(pfm) + 1), pfm, bg)

    rhits = rowMeans(as.matrix(rhits))
    return (list(fhits = as.vector(fhits), rhits = as.vector(rhits)))
}

#' Number of motif hits in a set of DNA sequences
#'
#' This function counts the number of motif hits that
#' are found in a given set of DNA sequences.
#'
#' Optionally, it can be used to count motif hits on
#' one or both strands, respectively.
#'
#'
#' @inheritParams scoreProfile
#' @param singlestranded Boolean that indicates whether a single strand or
#' both strands shall be scanned for motif hits.
#' Default: singlestranded = FALSE.
#' @return A list containing
#' \describe{
#' \item{nseq}{Number of individual sequences}
#' \item{lseq}{Vector of individual sequence lengths}
#' \item{numofhits}{Vector of the number of hits in each individual sequence}
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
#' # Count motif hits both strands
#' noc = motifcounter:::numMotifHits(seqs, motif, bg)
#' noc$numofhits
#'
#' # Count motif hits on a single strand
#' noc = motifcounter:::numMotifHits(seqs, motif, bg, singlestranded = TRUE)
#' noc$numofhits
#'
numMotifHits = function(seqs, pfm, bg, singlestranded = FALSE) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)

    stopifnot(class(seqs) == "DNAStringSet")

    # retrieve the number of motif hits
    noh = vapply(seqs, function(seq, pfm, bg, singlestranded) {
        ret = motifHits(seq, pfm, bg)
        if (singlestranded == FALSE) {
            return(sum(ret[[1]] + ret[[2]]))
        } else {
            return(sum(ret[[1]]))
        }
    }, integer(1), pfm, bg, singlestranded)

    # retrieve the individual sequence lengths
    # sequences containing "N" or "n" are assigned length zero
    lseq = lenSequences(seqs)

    nseq = length(seqs)

    return (list(
        nseq = nseq,
        lseq = lseq,
        numofhits = noh
    ))
}
