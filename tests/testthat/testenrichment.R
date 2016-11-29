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

    readBackground(seqfile,1)
    motif=t(as.matrix(read.table(motiffile)))

    motifEnrichmentTest(motif,seqs)
    motifEnrichmentTest(motif,seqs,singlestranded=TRUE)
    motifEnrichmentTest(motif,seqs,method="combinatorial")
    expect_error(motifEnrichmentTest(motif,seqs,singlestranded=TRUE,
                    method="combinatorial"))
    expect_error(motifEnrichmentTest(motif,seqs,singlestranded=TRUE,
                    method="tata"))
})
