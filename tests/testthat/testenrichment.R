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
            expect_equivalent(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)$numofhits,0)
            expect_equivalent(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)$pvalue,1)
            expect_true(is.na(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)$fold))
            expect_true(is.na(motifEnrichment(Biostrings::DNAStringSet(""),
                        motif,bg, method=method, singlestranded=ss)$logfold))
            expect_equivalent(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)$numofhits,0)
            expect_equivalent(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)$pvalue,1)
            expect_true(is.na(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)$fold))
            expect_true(is.na(motifEnrichment(
                        generateDNAStringSet(ncol(motif)-1,bg),
                        motif,bg, method=method, singlestranded=ss)$logfold))

            # check motif shorter than bg order
            bg=readBackground(seqs,2)
            expect_error(motifEnrichment(seqs,as.matrix(motif[, 1]),bg)$pvalue)
        }
    }

    # check wrong method
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="tata"))

    bg=readBackground(seqs,1)
    # run enrichment in different flavours
    m=motifEnrichment(seqs,motif,bg)
    expect_true(m$pvalue<1)
    expect_true(m$pvalue>=0)
    expect_true(is.numeric(m$fold))

    m=motifEnrichment(seqs,motif,bg,singlestranded=TRUE)
    expect_true(m$pvalue<1)
    expect_true(m$pvalue>=0)
    expect_true(is.numeric(m$fold))

    m=motifEnrichment(seqs,motif,bg,method="combinatorial")
    expect_true(m$pvalue<1)
    expect_true(m$pvalue>=0)
    expect_true(is.numeric(m$pvalue))

    # not supported
    expect_error(motifEnrichment(seqs,motif,bg,singlestranded=TRUE,
                    method="combinatorial"))

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs2=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs2,2)
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    
    
    expect_error(motifEnrichment(seqs,motif,bg))

    # DNAStringSet,motif
    expect_equivalent(motifEnrichment(seqs[1:2],motif,bg)$pvalue, 
                 0.1293926, tolerance=0.000001)
    
    # DNAStringSet,list
    ret <- motifEnrichment(seqs[1:2],
                 list(motif, motif),bg)$pvalue

    expect_equivalent(ret, matrix(rep(0.1293926, 2), 2, 1), tolerance=0.000001)

    ret <- motifEnrichment(Biostrings::DNAStringSetList(seqs[1:2],
                                                        seqs[1:2]),
                                 motif,bg)$pvalue
    # DNAStringSetList,matrix
    expect_equivalent(ret,
                 matrix(rep(0.1293926, 2), 1, 2), tolerance=0.000001)

    # DNAStringSetList,list
    expect_equal(motifEnrichment(Biostrings::DNAStringSetList(seqs[1:2],
                                                              seqs[1:2]),
                                 list(motif, motif),bg)$pvalue,
                 matrix(rep(0.1293926, 4), 2, 2, dimnames=list(NULL, NULL)), tolerance=0.000001)

})


test_that("motifEnrichment_w_ranges", {

    motifcounterOptions(alpha=0.001, gran=.1)

    # Load Sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    regions1 = system.file("extdata",
                           "example_jund.bed", package = "motifcounter")
    regions2 = system.file("extdata",
                           "example_arid3a.bed", package = "motifcounter")

    regions = GenomicRanges::GRangesList(
               "jund"=genomation::readBed(regions1),
               "arid3a"=genomation::readBed(regions2))
    
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    query = MotifDb::query
    
    motifs1 = as.list(query(query(MotifDb::MotifDb, "hsapiens"),
                      c("arid3a")))[1:3]
    motifs2 = as.list(query(query(MotifDb::MotifDb, "hsapiens"),
                      c("jund")))[1:3]
    motifs = c(motifs1, motifs2)
    
    # GRangesList,list
    result1 = motifEnrichment(regions, motifs, bg, genome=genome)
    
    expect_equal(dimnames(result1$pval), list(names(motifs), names(regions)))

    # list,list
    result2 = motifEnrichment(list(regions1, regions2), motifs, bg, genome=genome)
    expect_equivalent(result1$numofhits, result2$numofhits)
    expect_equivalent(result1$pvalue, result2$pvalue)
    expect_equivalent(result1$fold, result2$fold)
    expect_equivalent(result1$logfold, result2$logfold)

    # GRanges,list
    result3 = motifEnrichment(regions[[1]], motifs, bg, genome=genome)
    # character,list
    result4 = motifEnrichment(regions1, motifs, bg, genome=genome)

    expect_equivalent(result1$pvalue[,1], result3$pvalue)
    expect_equivalent(result1$pvalue[,1], result4$pvalue)

    # GRangesList,matrix
    result5 = motifEnrichment(regions, motifs[[1]], bg, genome=genome)
    expect_equivalent(result1$pvalue[1,], result5$pvalue)
    # list,matrix
    result6 = motifEnrichment(list(regions1, regions2), motifs[[1]], bg, genome=genome)
    expect_equivalent(result1$pvalue[1,], result6$pvalue)

    # GRanges,matrix
    result7 = motifEnrichment(regions[[1]], motifs[[1]], bg, genome=genome)
    # character,matrix
    result8 = motifEnrichment(regions1, motifs[[1]], bg, genome=genome)

    expect_equivalent(result1$pvalue[1,1], result7$pvalue)
    expect_equivalent(result1$pvalue[1,1], result8$pvalue)
    
})
