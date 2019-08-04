#' Generate DNAString
#'
#' This function generates a random DNAString of a given length
#' by sampling from the background model.
#'
#' @param len Integer length of the sequence
#' @template templates
#'
#' @return A DNAString object
#
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Generate a 1 kb random sequence
#' motifcounter:::generateDNAString(1000, bg)
#'
#' @seealso \code{\link{generateDNAStringSet}}
generateDNAString = function(len, bg) {
    len = as.integer(len)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    if (len < bg@order) {
        return(Biostrings::DNAString(""))
    }
    seq = paste(rep("N", len), collapse = "", sep = "")
    ret = .C(
        motifcounter_generateRndSeq,
        as.character(seq),
        as.integer(len),
        getStation(bg),
        getTrans(bg),
        as.integer(getOrder(bg))
    )
    return(Biostrings::DNAString(ret[[1]]))
}

#' Generate DNAStringSet
#'
#' This function generates a DNAStringSet-object of the
#' given individual sequence lengths
#' by sampling from the background model.
#'
#'
#' @inheritParams compoundPoissonDist
#' @template templates
#'
#' @return A DNAStringSet object
#
#' @examples
#'
#' # Load sequences
#' seqfile = system.file("extdata", "seq.fasta", package = "motifcounter")
#' seqs = Biostrings::readDNAStringSet(seqfile)
#'
#' # Load background
#' bg = readBackground(seqs, 1)
#'
#' # Generate random sequences of various lengths
#' motifcounter:::generateDNAStringSet(10:50, bg)
#'
#' @seealso \code{\link{generateDNAStringSet}}
generateDNAStringSet = function(seqlen, bg) {
    stopifnot(is(bg, "Background"))
    validObject(bg)
    seqs = c()
    for (i in seq_len(length(seqlen))) {
        seqs = c(seqs, generateDNAString(seqlen[i], bg))
    }

    return(Biostrings::DNAStringSet(seqs))
}

#' Empirical number of motif hits distribution
#'
#' This function repeatedly simulates random DNA sequences according to
#' the background model and
#' subsequently counts how many motif hits occur in them.
#' Thus, this function gives rise to the empirical
#' distribution of the number of motif hits.
#' This function is only used for benchmarking analysis.
#'
#' @inheritParams numMotifHits
#' @inheritParams compoundPoissonDist
#' @param nsim Integer number of random samples.
#' @return A List that contains
#' \describe{
#' \item{dist}{Empirical distribution of the number of motif hits}
#' }
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
#' # Study the counts in one sequence of length 150 bp
#' seqlen = rep(150, 1)
#'
#' # Compute empirical distribution of the number of motif hits
#' # by scanning both strands using 100 samples
#' simc = motifcounter:::simulateNumHitsDist(motif, bg,
#'     seqlen, nsim = 100, singlestranded = FALSE)
#'
#' # Compute empirical distribution of the number of motif hits
#' # by scanning a single strand using 100 samples
#' simc = motifcounter:::simulateNumHitsDist(motif, bg,
#'     seqlen, nsim = 100, singlestranded = TRUE)
#'
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
simulateNumHitsDist = function(pfm, bg, seqlen, nsim, singlestranded = FALSE) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)
    stopifnot(nsim > 0)

    freq = rep(0, 10)
    for (i in seq_len(nsim)) {
        seqs = generateDNAStringSet(seqlen, bg)
        nom = numMotifHits(seqs, pfm, bg, singlestranded)
        nom = sum(nom$numofhits)
        if (nom >= length(freq)) {
            #expand freq
            new_freq = rep(0, nom + 1)
            new_freq[seq_len(length(freq))] = freq
            freq = new_freq
        }
        freq[nom + 1] = freq[nom + 1] + 1
    }
    freq = freq / sum(freq)
    return(list(dist = freq))
}

#' Empirical clump size distribution
#'
#' This function repeatedly simulates random DNA sequences according to
#' the background model and
#' subsequently counts the number of k-clump occurrences, where
#' denotes the clump size.
#' This function is only used for benchmarking analysis.
#'
#' @inheritParams numMotifHits
#' @inheritParams compoundPoissonDist
#' @param nsim Integer number of random samples.
#' @return A List that contains
#' \describe{
#' \item{dist}{Empirical distribution of the clump sizes}
#' }
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
#' # Study the clump size frequencies in one sequence of length 1 Mb
#' seqlen = 1000000
#'
#' # scan both strands
#' simc = motifcounter:::simulateClumpSizeDist(motif, bg, seqlen)
#'
#' # scan a single strand
#' simc = motifcounter:::simulateClumpSizeDist(motif, bg,
#'     seqlen, singlestranded = TRUE)
#'
#' @seealso \code{\link{compoundPoissonDist}},\code{\link{combinatorialDist}}
simulateClumpSizeDist = function(pfm, bg, seqlen,
                                nsim=10, singlestranded = FALSE) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)
    stopifnot (length(seqlen) == 1)

    # 
    freq = rep(0, 5)
    for (p in seq_len(nsim)) {
        seq = generateDNAString(seqlen, bg)
        hits = motifHits(seq, pfm, bg)

        # determine the clump sizes
        if (singlestranded == TRUE) {
            h=hits$fhits
        } else {
            h=hits$fhits + hits$rhits
        }

        index_sequence = which(h>0)
        previous_index = index_sequence[1]
        cnt=h[previous_index]
        for (i in tail(seq_len(length(index_sequence)),-1)) {
            if (index_sequence[i] < previous_index+ncol(pfm)) {
                cnt = cnt + h[index_sequence[i]]
            } else {
                if (length(freq) < cnt) {
                    x = rep(0, cnt)
                    x[seq_len(length(freq))] = freq
                    freq = x
                }
                freq[cnt] = freq[cnt] + 1
                cnt = h[index_sequence[i]]
            }
            previous_index = index_sequence[i]
        }
    }

    freq = freq / sum(freq)
    return(list(dist = freq))
}

#' Empirical score distribution
#'
#' This function estimates the empirical score distribution
#' on a set of randomly generated DNA sequences based on the
#' background model.
#' This function is only used for benchmarking analysis.
#'
#'
#' @inheritParams simulateNumHitsDist
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
#'
#' # Compoute the empirical score distribution in
#' # sequences of length 1kb using 1000 samples
#' motifcounter:::scoreDistEmpirical(motif, bg, seqlen = 1000, nsim = 1000)
#'
#' @seealso \code{\link{scoreDist}}
scoreDistEmpirical = function(pfm, bg, seqlen, nsim) {
    motifValid(pfm)
    stopifnot(is(bg, "Background"))
    validObject(bg)
    motifAndBackgroundValid(pfm, bg)
    stopifnot(nsim > 0)
    if (seqlen < ncol(pfm)) {
        warning("seqlen shorter than ncol(pfm). Scores will always be zero.")
    }

    seqs = generateDNAStringSet(rep(seqlen, nsim), bg)
    sh = scoreHistogram(seqs, pfm, bg)

    probs = sh$dist / sum(sh$dist)
    return(list(scores = sh$score, dist = probs))
}
