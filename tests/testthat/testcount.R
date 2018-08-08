context("Count")
library(MotifDb)

test_that("oct4_vignette_example_alpha0001", {
  library(MotifDb)
  alpha=0.001
  gran=0.1
  motifcounterOptions(alpha, gran)

  order <- 1
  file <- system.file("extdata", "seq.fasta", package = "motifcounter")
  seqs <- Biostrings::readDNAStringSet(file)
  bg <- readBackground(seqs, order)

  file <- system.file("extdata", "oct4_chipseq.fa", package = "motifcounter")
  oct4peaks <- Biostrings::readDNAStringSet(file)

  oct4 <- as.list(query(query(query(MotifDb, "hsapiens"),
                              "pou5f1"), "jolma2013"))[[1]]
  motif <- oct4

  mhits <- motifHits(oct4peaks[[1]], motif, bg)

  # Inspect the result
  fhitpos <- which(mhits$fhits == 1)
  rhitpos <- which(mhits$rhits == 1)
  expect_equal(length(fhitpos), 0)
  expect_equal(length(rhitpos), 1)
  expect_equal(rhitpos, 94)

})

test_that("oct4_vignette_example_alpha001", {
  library(MotifDb)
  alpha=0.01
  gran=0.1
  motifcounterOptions(alpha, gran)

  order <- 1
  file <- system.file("extdata", "seq.fasta", package = "motifcounter")
  seqs <- Biostrings::readDNAStringSet(file)
  bg <- readBackground(seqs, order)

  file <- system.file("extdata", "oct4_chipseq.fa", package = "motifcounter")
  oct4peaks <- Biostrings::readDNAStringSet(file)

  oct4 <- as.list(query(query(query(MotifDb, "hsapiens"),
                              "pou5f1"), "jolma2013"))[[1]]
  motif <- oct4

  mhits <- motifHits(oct4peaks[[1]], motif, bg)

  # Inspect the result
  fhitpos <- which(mhits$fhits == 1)
  rhitpos <- which(mhits$rhits == 1)
  expect_equal(length(fhitpos), 4)
  expect_equal(length(rhitpos), 4)
  expect_equal(fhitpos, c(54, 87, 93, 112))
  expect_equal(rhitpos, c(55, 94, 111, 118))
  
})

test_that("motifHits", {
    # init 
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # load sequences
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    # Load motif
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))


    # Check too short sequence
    expect_equal(length(motifHits(Biostrings::DNAString(""),motif,bg)[[1]]),0)
    expect_equal(length(motifHits(Biostrings::DNAString(""),motif,bg)[[2]]),0)
    expect_equal(length(motifHits(generateDNAString(ncol(motif)-1,bg),motif,bg)[[1]]),0)
    expect_equal(length(motifHits(generateDNAString(ncol(motif)-1,bg),motif,bg)[[2]]),0)

    # check motif shorter than bg order
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,2)
    expect_error(motifHits(seqs,as.matrix(motif[, 1]),bg)[[1]])

    # Reload background
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    # only use the first two sequences without N's to retain the current testing behaviour
    bg=readBackground(seqs[1:2],1)

    # Test correct hit positions
    mh=motifHits(seqs[[1]], motif,bg)
    # because it is a palindrome
    expect_equal(mh[[1]],mh[[2]])
    # there must be a hit at position 13
    x=rep(0,20)
    x[13]=1  
    expect_equal(mh[[1]],x)
    
    # no hit in the second sequence
    mh=motifHits(seqs[[2]], motif,bg)
    expect_equal(mh[[1]],mh[[2]])
    x=rep(0,length(seqs[[2]])-ncol(motif)+1)
    expect_equal(mh[[1]],x)

    # zero length sequence, due to 'N'
    mh=motifHits(seqs[[3]], motif,bg)
    expect_equal(c(length(mh[[1]]),length(mh[[2]])),c(10,10))
    expect_equal(mh$fhits, c(0, NaN, NaN, NaN, NaN, NaN, NaN, 0, 0, 0))

    # Check with aaa repeat motif
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    mh=motifHits(seqs[[2]],motif,bg)
    x=rep(0,length(seqs[[2]])-ncol(motif)+1)
    # in the reverse strand no hit
    expect_equal(mh[[2]],x)
    x[4]=1
    # in the forward strand one hit
    expect_equal(mh[[1]],x)

})

test_that("motifHitProfile", {
    # init
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Load seqs
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load Background
    bg=readBackground(seqs[1:2],1)
    
    # palindrom
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    # error because it is no DNAStringSet
    expect_error(motifHitProfile(seqs[[1]],motif,bg)) 

    # error because unequal sequence length
    expect_error(motifHitProfile(seqs[1:2],motif,bg))

    # Check too short sequence
    expect_equal(length(motifHitProfile(Biostrings::DNAStringSet(""),
                                        motif,bg)$fhits),0)
    expect_equal(length(motifHitProfile(generateDNAStringSet(ncol(motif)-1,bg),
                                 motif,bg)$fhits),0)

    # check motif shorter than bg order
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs2=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs2,2)
    expect_error(motifHitProfile(seqs,as.matrix(motif[, 1]),bg)[[1]])

    # Reload background
    bg=readBackground(seqs[1:2],1)

    # Test correct hit positions
    mh=motifHitProfile(seqs[1], motif,bg)
    # because it is a palindrome
    expect_equal(mh[[1]],mh[[2]])
    # there must be a hit at position 13
    x=rep(0,20)
    x[13]=1  
    expect_equal(mh[[1]],x)
    
    # no hit in the second sequence
    mh=motifHitProfile(seqs[2], motif,bg)
    expect_equal(mh[[1]],mh[[2]])
    x=rep(0,length(seqs[[2]])-ncol(motif)+1)
    expect_equal(mh[[1]],x)

    # zero length sequence, due to 'N'
    mh=motifHitProfile(seqs[3], motif,bg)
    expect_equal(c(length(mh[[1]]),length(mh[[2]])),c(10,10))
    expect_equal(mh$fhits, c(0, NaN, NaN, NaN, NaN, NaN, NaN, 0, 0, 0))
    
    seqs2=seqs
    seqs2[[1]] = seqs2[[1]][1:10]
    seqs2[[3]] = seqs2[[3]][1:10]
    mh=motifHitProfile(seqs2, motif,bg)
    
    # Check with aaa repeat motif
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))


    profile=motifHitProfile(seqs[2],motif,bg)
    x=rep(0,length(seqs[[2]])-ncol(motif)+1)
    # in the reverse strand no hit
    expect_equal(profile[[2]],x)
    x[4]=1
    # in the forward strand one hit
    expect_equal(profile[[1]],x)

})

test_that("numMotifHits", {
    # init
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Load seqs
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load Background
    bg=readBackground(seqs[1:2],1)
    
    # palindrom
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))


    # no DNAString
    expect_error(numMotifHits(seqfile,motif,bg,singlestranded=TRUE))
    expect_error(numMotifHits(seq[[1]],motif,bg,singlestranded=TRUE))

    # Check short DNA sequence
    for (ss in c(TRUE, FALSE)) {

        expect_equal(numMotifHits(Biostrings::DNAStringSet(""),motif,bg)[[1]],1)
        expect_equal(numMotifHits(Biostrings::DNAStringSet(""),motif,bg)[[2]],0)
        expect_equal(numMotifHits(Biostrings::DNAStringSet(""),motif,bg)[[3]],0)
        expect_equal(numMotifHits(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)[[1]],1)
        expect_equal(numMotifHits(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)[[2]],ncol(motif)-1)
        expect_equal(numMotifHits(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)[[3]],0)
    }


    
    # test with palindromic motif
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    # there is a hit in the first sequence
    nom=numMotifHits(seqs,motif,bg,singlestranded=TRUE)
    expect_equal(nom$nseq,3)
    expect_equal(as.vector(nom$lseq),c(23,10,13))
    expect_equal(as.vector(nom$numofhits),c(1,0,NaN))

    nom=numMotifHits(seqs,motif,bg,singlestranded=FALSE)
    expect_equal(nom$nseq,3)
    expect_equal(as.vector(nom$lseq),c(23,10,13))
    expect_equal(as.vector(nom$numofhits),c(2,0,NaN))


    # test with repeat motif
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    # there is a hit in the second sequence
    nom=numMotifHits(seqs,motif,bg,singlestranded=TRUE)
    expect_equal(nom$nseq,3)
    expect_equal(as.vector(nom$lseq),c(23,10,13))
    expect_equal(as.vector(nom$numofhits),c(0,1,NaN))

    nom=numMotifHits(seqs,motif,bg,singlestranded=FALSE)
    expect_equal(nom$nseq,3)
    expect_equal(as.vector(nom$lseq),c(23,10,13))
    expect_equal(as.vector(nom$numofhits),c(0,1,NaN))

})

