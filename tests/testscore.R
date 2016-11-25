
library(motifcounter)
alpha=0.01

motifcounterOption(alpha)
pwmname="x1.tab"
pwmfile=system.file("extdata",pwmname, package="motifcounter")
seqfile=system.file("extdata","seq.fasta", package="motifcounter")

# 1. test:
# with an order zero model only sample one letter for P times
# and compare to the actual probability

readMotif(pwmfile)
pwm=motif2matrix()
readMotif(as.matrix(pwm[,1]),0)
pwm=motif2matrix()
seqlen=ncol(pwm)

readBackground(seqfile, 0)
readBackgroundForSampling(seqfile,0)
# compute the score distribution for this case in R
bg=fetchBackground()
s=round((log(pwm)-log(bg[[1]]))*10)
srange=seq(min(s),max(s))
p=rep(0,length(srange))
for (i in 1:length(srange)){
    p[i]=sum(bg[[2]][which(s==srange[i])])
}

dp=scoreDist()

if (!all(srange==round(dp[[1]]*10))) {
    stop("Score ranges between dynamic programming 
        variant and R variant divergent")
}
if (!all(abs(p-dp[[2]])<1e-12)) {
    stop("Score probabilities between dynamic 
        programming variant and R variant divergent")
}

# test whether the range is equally long
# test whether the zero entries in the score distribution overlap perfectly
# test with the stationary probabilities of the background model only

for (m in seq(0,3)) {
    print(m)
    readBackground(seqfile, m)
    readBackgroundForSampling(seqfile,m)
    readMotif(pwmfile)
    pwm=motif2matrix()
    if (m==0) {
        readMotif(as.matrix(pwm[,1]),0)
    } else {
        readMotif(as.matrix(pwm[,1:m]),0)
    }
    pwm=motif2matrix()
    seqlen=ncol(pwm)

    #simluate score distribution
    sims=simulateScoreDist(seqlen,100000)
    #compute exact score distribution
    dp=scoreDist()
    bf=scoreDistBf()

    if (!all(bf[[1]]==dp[[1]])) {
      stop(paste("Score ranges between dynamic programming 
           variant and enumerative variant divergent: ",m))
    }
    if (!all(bf[[2]]==dp[[2]])) {
      stop(paste("Score distribution between dynamic programming 
           variant and enumerative variant divergent: ",m))
    }
    if (!all(sims[[1]]==dp[[1]])) {
      stop(paste("Score ranges between dynamic programming 
           variant and sampling variant divergent: ", m))
    }
    if (!all(!xor(sims[[2]]==0,dp[[2]]==0))) {
      stop(paste("Score probabilities between dynamic 
           programming variant and simulated variant divergent: ",m))
    }
    # This test is incorrect
    # Due to rounding differences of scores collected on either
    # DNA strand, this condition might actually be wrong
    #if (dp[[2]][1]<=0 || dp[[2]][length(dp[[2]])]<=0) {
    #  stop(paste("The first and the last 
    #             score entry must be greater than zero: ",m))
    #}
}

# test whether the range is equally long
# test whether the zero entries in the score distribution overlap perfectly
# test with the stationary and the transition probabilities
for (m in seq(1,3)) {
    print(m)
    readBackground(seqfile, m)
    readBackgroundForSampling(seqfile,m)
    readMotif(pwmfile)
    pwm=motif2matrix()
    readMotif(as.matrix(pwm[,1:(m+1)]),0)
    pwm=motif2matrix()
    seqlen=ncol(pwm)

    #simluate score distribution
    sims=simulateScoreDist(seqlen,1000000)
    #compute exact score distribution
    dp=scoreDist()
    bf=scoreDistBf()

    if (!all(bf[[1]]==dp[[1]])) {
      stop(paste("Score ranges between dynamic programming 
           variant and enumerative variant divergent: ",m))
    }
    if (!all(abs(bf[[2]]-dp[[2]])<1e-5)) {
      stop(paste("Score distribution between dynamic programming 
           variant and enumerative variant divergent: ",m))
    }
    if (!all(sims[[1]]==dp[[1]])) {
      stop(paste("Score ranges between dynamic programming 
           variant and sampling variant divergent: ", m))
    }
    if (!all(!xor(sims[[2]]==0,dp[[2]]==0))) {
      stop(paste("Score probabilities between dynamic 
           programming variant and simulated variant divergent: ",m))
    }
    # This test is incorrect
    # Due to rounding differences of scores collected on either
    # DNA strand, this condition might actually be wrong
    #if (dp[[2]][1]<=0 || dp[[2]][length(dp[[2]])]<=0) {
      #stop(paste("The first and the last 
                 #score entry must be greater than zero: ",m))
    #}
}

