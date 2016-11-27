context("Combinatorial model")

test_that("combinatorial", {
    alpha=0.01
    gran=0.1
    seqlen=100
    numofseqs=10
    maxhits=100
    motifcounterOption(alpha, gran)

    pwmname="x3.tab"
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    motiffile=system.file("extdata",pwmname, package="motifcounter")

    readBackground(seqfile,1)
    motif=t(as.matrix(read.table(motiffile)))

    op=probOverlapHit(motif)

    dist=combinatorialDist(seqlen, op)

    plot(dist$dist)
})
