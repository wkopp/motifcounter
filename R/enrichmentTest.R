#' Internal logic Enrichment of motif hits
#'
#' This function determines whether a given motif is enriched in a given
#' DNA sequences.
#'
#' Enrichment is tested by comparing the observed
#' number of motif hits against a theoretical distribution of the number
#' of motif hits in random DNA sequences.
#' Optionally, the theoretical distribution of the number of motif
#' hits can be evaluated by either a 'compound Poisson model'
#' or the 'combinatorial model'.
#' Additionally, the enrichment test can be conducted with respect
#' to scanning only the forward strand or both strands of the DNA
#' sequences. The latter option is only available for the
#' 'compound Poisson model'.
#' @include count_wrapper.R
#'
#' @inheritParams numMotifHits
#' @param method String that defines whether to use
#' the 'compound' Poisson approximation' or the 'combinatorial' model.
#' Default: method='compound'.
#' @param ... externally supplied overlap object
#'
#' @return List that contains
#' \describe{
#' \item{numofhits}{Number of observations}
#' \item{pvalue}{P-value for the enrichment test}
#' \item{fold}{Fold-enrichment with respect to the expected number of hits}
#' \item{logfold}{log-Fold-enrichment with respect to
#'                the expected number of hits}
#' }
#' @seealso \code{\link{compoundPoissonDist}}, \code{\link{combinatorialDist}}
motifEnrichment_ = function(seqs, pfm, bg, singlestranded = FALSE,
                            method = "compound", ...) {

    #print("motifEnrichment_")
    stopifnot(is(bg, "Background"))
    motifValid(pfm)

    validObject(bg)
    motifAndBackgroundValid(pfm, bg)
    stopifnot(is.logical(singlestranded))

    #compute overlapping hit probs
    if ("overlap" %in% names(list(...))) {
        overlap <- list(...)$overlap
    } else {
        overlap = probOverlapHit(pfm, bg, singlestranded)
    }

    # detemine the number of motif hits
    observations = numMotifHits(seqs, pfm, bg, singlestranded)

    if (method == "compound") {
        dist = compoundPoissonDist(observations$lseq, overlap)
    } else if (method == "combinatorial") {
        dist = combinatorialDist(observations$lseq, overlap)
    } else {
        stop("Invalid method: 'method' must be 'compound' or 'combinatorial'")
    }
    ind = seq_len(length(dist$dist))
    ind = ind[ind >= (sum(observations$numofhits) + 1)]
    p = sum(dist$dist[ind])

    return (list(
        numofhits = sum(observations$numofhits),
        pvalue = p,
        fold = (sum(observations$numofhits) +1)/
            (sum(dist$dist * seq(0, length(dist$dist) - 1)) + 1),
        logfold = log((sum(observations$numofhits) + 1) /
            (sum(dist$dist * seq(0, length(dist$dist) - 1))+1)),
        meanhits = sum(observations$numofhits) / sum(observations$lseq)
    ))
}

#' Enrichment of motif hits
#'
#' This function determines whether a given motif is enriched in a given
#' DNA sequences.
#'
#' Enrichment is tested by comparing the observed
#' number of motif hits against a theoretical distribution of the number
#' of motif hits in random DNA sequences.
#' Optionally, the theoretical distribution of the number of motif
#' hits can be evaluated by either a 'compound Poisson model'
#' or the 'combinatorial model'.
#' Additionally, the enrichment test can be conducted with respect
#' to scanning only the forward strand or both strands of the DNA
#' sequences. The latter option is only available for the
#' 'compound Poisson model'.
#' @name motifEnrichment
#' @aliases motifEnrichment
#' @param seqs_or_regions DNAStringSet, DNAStringSetList,
#'                        GRanges, GRangesList,
#'                        character-string denoting bed-file location or
#'                        list(character) denoting a list of bed-files
#' @param pfms matrix or list of matrices representing position weight matrices
#' @param normalize_pfm If TRUE, normalization of the PFM is performed
#'                        to avoid zero-entries. Default: FALSE
#' @param genome BSgenome object from which the sequences are extracted
#'               using the provided regions. Only used in conjunction with
#'               GRanges, GRangesList, and regions extracted from bed-files.
#' @include count_wrapper.R
#'
#' @inheritParams numMotifHits
#' @param method String that defines whether to use
#' the 'compound' Poisson approximation' or the 'combinatorial' model.
#' Default: method='compound'.
#' @param ... see concrete methods
#'
#' @return List that contains
#' \describe{
#' \item{pvalue}{P-value for the enrichment test.}
#' \item{fold}{Fold-enrichment with respect to the expected number of hits}
#' \item{logfold}{log-Fold-enrichment with respect to
#'                the expected number of hits}
#' }
#' @import BSgenome
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
#' # 1 ) Motif enrichment test w.r.t. scanning a *single* DNA strand
#' # based on the 'Compound Poisson model'
#'
#' result = motifEnrichment(seqs, motif, bg,
#'             singlestranded = TRUE, method = "compound")
#'
#' # 2 ) Motif enrichment test w.r.t. scanning *both* DNA strand
#' # based on the 'Compound Poisson model'
#'
#' result = motifEnrichment(seqs, motif, bg, method = "compound")
#'
#' # 3 ) Motif enrichment test w.r.t. scanning *both* DNA strand
#' # based on the *combinatorial model*
#'
#' result = motifEnrichment(seqs, motif, bg, singlestranded = FALSE,
#'             method = "combinatorial")
#'
#' # 4 ) Evaluate the enrichment in two sets of sequences
#' #     (e.g. promoters and enhancers) using three motifs.
#'
#' result = motifEnrichment(Biostrings::DNAStringSetList(seqs, seqs),
#'                          list(motif, motif, motif),
#'                          bg)
#'
#' # 4 ) Evaluate the enrichment in two sets of genomic regions
#' #     exracted from the human genome
#' #     (e.g. promoters and enhancers) using three motifs
#'
#' regions1 = system.file("extdata",
#'                        "example_jund.bed", package = "motifcounter")
#' regions2 = system.file("extdata",
#'                        "example_arid3a.bed", package = "motifcounter")
#' regions = GenomicRanges::GRangesList(
#'            "jund"=genomation::readBed(regions1),
#'            "arid3a"=genomation::readBed(regions2))
#'
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' query <- MotifDb::query
#' motifs1 = as.list(query(query(MotifDb::MotifDb, "hsapiens"), c("arid3a")))
#' motifs2 = as.list(query(query(MotifDb::MotifDb, "hsapiens"), c("jund")))
#' motifs = c(motifs1[1], motifs2[1])
#'
#' result = motifEnrichment(regions, motifs, bg, genome=genome)
#'
#' @seealso \code{\link{compoundPoissonDist}}, \code{\link{combinatorialDist}}
setGeneric(
    name="motifEnrichment",
    def = function(seqs_or_regions, pfms, ...) {
        standardGeneric("motifEnrichment")
    }
)

#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
setMethod("motifEnrichment",
    signature(seqs_or_regions="DNAStringSet", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome = NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("DNAStringSet,matrix")
        stopifnot(is.null(genome))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        # with motifEnrichment_: returns list(scalar)
        return(motifEnrichment(DNAStringSetList(seqs_or_regions), pfms, bg,
            singlestranded=singlestranded,
            method=method,
            normalize_pfm=normalize_pfm))
    }
)

#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
setMethod("motifEnrichment",
    signature(seqs_or_regions="DNAStringSetList", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("DNAStringSetList,matrix")
        stopifnot(is.null(genome))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        if (normalize_pfm) {
            pfms = normalizeMotif(pfms)
        }

        overlap = probOverlapHit(pfms, bg, singlestranded)

        res = lapply(seqs_or_regions, function(seqs, pfms, bg, singlestranded,
                    method,
                    normalize_pfm, overlap)
                motifEnrichment_(seqs, pfms, bg, genome=NULL,
                    singlestranded=singlestranded,
                    method=method,
                    overlap=overlap),
                pfms, bg, singlestranded, method, normalize_pfm,
                overlap)

        # with DNAStringSet,matrix: list(vector)
        return(collapseListRegions(res, length(res), names(seqs_or_regions)))

    }
)



#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
#' @importFrom GenomicRanges GRanges
setMethod("motifEnrichment",
    signature(seqs_or_regions="GRanges", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("GRanges,matrix")

        stopifnot(is(genome, "BSgenome"))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        #seqs = getSeq(genome, seqs_or_regions)

        # with DNAStringSet,matrix: returns list(scalar values)
        return(motifEnrichment(GRangesList(seqs_or_regions), pfms, bg,
                            genome=genome,
                            singlestranded=singlestranded,
                            method=method,
                            normalize_pfm=normalize_pfm))
    }
)


#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
#' @importFrom GenomicRanges GRangesList
setMethod("motifEnrichment",
    signature(seqs_or_regions="GRangesList", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {

        #print("GRangesList,matrix")
        stopifnot(is(genome, "BSgenome"))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        seqslist = DNAStringSetList(lapply(seqs_or_regions,
            function(x, genome) {
                return(getSeq(genome, x))
            }, genome))
        names(seqslist) <- names(seqs_or_regions)

        res = motifEnrichment(seqslist, pfms, bg, genome=NULL,
                singlestranded=singlestranded,
                method=method,
                normalize_pfm=normalize_pfm)
        # with DNAStringSetList,matrix: list(vector)
        return(collapseListRegions(res, length(res), names(seqslist)))

    }
)


#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
#' @importFrom genomation readBed
setMethod("motifEnrichment",
    signature(seqs_or_regions="character", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("character,matrix")
        stopifnot(is(genome, "BSgenome"))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        # with list,matrix: returns list(vector)
        res = motifEnrichment(list(seqs_or_regions), pfms, bg,
            genome=genome, singlestranded=singlestranded,
            method=method,
            normalize_pfm=normalize_pfm)
        return(collapseListRegions(res, length(res), NULL))
    }
)

#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
#' @importFrom genomation readBed
setMethod("motifEnrichment",
    signature(seqs_or_regions="list", pfms="matrix"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("list,matrix")

        stopifnot(is(genome, "BSgenome"))
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))
        stopifnot(all(file.exists(unlist(seqs_or_regions))))

        regions = GRangesList(lapply(seqs_or_regions, readBed))
        names(regions) = basename(unlist(seqs_or_regions))

        # with GRangesList,matrix: returns list(vector)
        res = motifEnrichment(regions, pfms, bg, genome=genome,
            singlestranded=singlestranded,
            method=method,
            normalize_pfm=normalize_pfm)
        return(collapseListRegions(res, length(res), names(regions)))
    }
)


#' @rdname motifEnrichment
#' @aliases motifEnrichment
#' @export
setMethod("motifEnrichment",
    signature(pfms="list"),
    function(seqs_or_regions, pfms, bg, genome=NULL, singlestranded = FALSE,
            method = "compound",
            normalize_pfm=TRUE) {
        #print("*,list")
        stopifnot(is.logical(singlestranded))
        stopifnot(is.logical(normalize_pfm))

        res = lapply(pfms, function(motif, seqs_or_regions, bg, genome,
                                    singlestranded,
                                    method, normalize_pfm)
                                    motifEnrichment(seqs_or_regions, motif, bg,
                                        genome=genome,
                                        singlestranded=singlestranded,
                                        method=method,
                                        normalize_pfm=normalize_pfm),
                                seqs_or_regions, bg, genome, singlestranded,
                                method, normalize_pfm)

        # with any,matrix: returns list(matrix)
        return(collapseListMotifs(res, length(pfms), names(pfms)))

    }
)

collapseListRegions <- function(res, nregions, regionnames) {
    # this function collates results in a list(list(scalar))
    # to list(vector).
    # the result summarizes the results for a single motif over
    # multiple regions.

    if (is(res, "list") && length(intersect(c("numofhits",
        "pvalue", "fold", "logfold", "meanhits"), names(res)))==5) {
        return (res)
    }

    dnames <- list(NULL, regionnames)

    newres = list(numofhits=matrix(0, nrow=1, ncol=nregions, dimnames=dnames),
                pvalue=matrix(0, nrow=1, ncol=nregions, dimnames=dnames),
                fold=matrix(0, nrow=1, ncol=nregions, dimnames=dnames),
                logfold=matrix(0, nrow=1, ncol=nregions, dimnames=dnames),
                meanhits=matrix(0, nrow=1, ncol=nregions, dimnames=dnames))

    # regions are always first collected in the inner loop
    # that is from the scalar values we obtain a vector
    for (k in seq_len(length(newres))) {
        for (i in seq_len(nregions)) {
            newres[[k]][1, i] = res[[i]][[k]]
        }
    }

    return(newres)
}

collapseListMotifs <- function(res, nmotifs, motifnames) {
    # this function collates results in a list(list(vector))
    # to list(matrix).
    # it is applied to the results of multiple motifs on
    # multple sequences.

    #print("collapseListMotifs")
    if (is(res, "list") && length(intersect(c("numofhits",
        "pvalue", "fold", "logfold", "meanhits"), names(res)))==5) {
        return (res)
    }

    # extract the number of regions and region names
    # from the results so far.
    nregions <- length(res[[1]][[2]])
    dnames <- list(motifnames, dimnames(res[[1]][[1]])[[2]])

    newres = list(numofhits=matrix(0, nrow=nmotifs,
            ncol=nregions, dimnames=dnames),
        pvalue=matrix(0, nrow=nmotifs, ncol=nregions, dimnames=dnames),
        fold=matrix(0, nrow=nmotifs, ncol=nregions, dimnames=dnames),
        logfold=matrix(0, nrow=nmotifs, ncol=nregions, dimnames=dnames),
        meanhits=matrix(0, nrow=nmotifs, ncol=nregions, dimnames=dnames))

    for (k in seq_len(length(newres))) {
        for (j in seq_len(nmotifs)) {
            newres[[k]][j,] = res[[j]][[k]]
        }
    }

    return(newres)
}
