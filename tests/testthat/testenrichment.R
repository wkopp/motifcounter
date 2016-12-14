context("Enrichment")

test_that("motifEnrichment", {
    # init
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Load Sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    # Load motif
    pwmname="x3.tab"
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # test short sequence
    for (method in c("compound", "combinatorial")) {
         for (ss in c(TRUE,FALSE)) {
            bg=readBackground(seqs,1)
            expect_equal(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)[[1]],1)
            expect_true(is.na(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)[[2]]))
            expect_equal(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)[[1]],1)
            expect_true(is.na(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)[[2]]))

            # check motif shorter than bg order
            bg=readBackground(seqs,2)
            expect_error(motifEnrichment(seqs,as.matrix(motif[, 1]),bg)[[1]])
        }
    }

    # check wrong method
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="tata"))

    bg=readBackground(seqs,1)
    # run enrichment in different flavours
    m=motifEnrichment(seqs,motif,bg)
    expect_true(m[[1]]<1)
    expect_true(m[[1]]>=0)
    expect_true(is.numeric(m[[2]]))

    m=motifEnrichment(seqs,motif,bg,singlestranded=TRUE)
    expect_true(m[[1]]<1)
    expect_true(m[[1]]>=0)
    expect_true(is.numeric(m[[2]]))

    m=motifEnrichment(seqs,motif,bg,method="combinatorial")
    expect_true(m[[1]]<1)
    expect_true(m[[1]]>=0)
    expect_true(is.numeric(m[[2]]))

    # not supported
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="combinatorial"))
})
