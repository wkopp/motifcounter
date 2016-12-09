context("Enrichment test")

test_that("enrichment", {
    alpha=0.01
    gran=0.1
    seqlen=100
    numofseqs=10
    motifcounterOption(alpha, gran)

    pwmname="x3.tab"
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    motiffile=system.file("extdata",pwmname, package="motifcounter")

    seqs=Biostrings::readDNAStringSet(seqfile)

    bg=readBackground(seqs,1)
    motif=t(as.matrix(read.table(motiffile)))

    motifEnrichment(seqs,motif,bg)
    motifEnrichment(seqs,motif,bg,singlestranded=TRUE)
    motifEnrichment(seqs,motif,bg,method="combinatorial")
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="combinatorial"))
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="tata"))
})
