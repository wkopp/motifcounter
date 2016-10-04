library(devtools)
library(mdist)
library(seqLogo)
library(ggplot2)
paperdir=Sys.getenv("CPPAPER")
paperdir="/project/Tf2motif/home/thesis"
figdir=sprintf("%s/figures", paperdir)
variants=data.frame(
                markov=c(rep(1,1),rep(2,1)),
                name=rep(c("sp1"),2),
                type=rep(c("transfac"),2))

scorelabnames1=c(expression("P"["E,1"]*"(S)"), expression("P"["E,2"]*"(S)"))
scorelabnames2=c(expression("P"["1"]*"(S)"), expression("P"["2"]*"(S)"))
countlabnames1=c(expression("P"["E,1"]*"(X)"),expression("P"["E,2"]*"(X)"))
countlabnames2=c(expression("P"["CP,1"]^"N"*"(X)"), expression("P"["CP,2"]^"N"*"(X)"))
alpha=0.01
lalpha=-log10(alpha)
gran=0.1
mdist.option(alpha, gran)
seqlen=10000
numofseqs=1
maxhits=500

seqfile=system.file("extdata","cpg.fa", package="mdist")
#print.background()
#print.background.sampling()
#reldiff=rep(0,nrow(variants))
#rel2=reldiff
reldiff=rep(0,nrow(variants))
#reldiff[1:12]=rel2;

for (var in 1:nrow(variants)) {
    markov=variants$markov[var]
    read.background(seqfile,markov)
    pwmname=variants$name[var]
    read.background.sampling(seqfile,variants$markov[var])
    print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,
                                                           variants$markov[var]))

    pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
    package="mdist")

    #load motif model from file
    read.motif(pwmfile,variants$type[var], 0.01)
    pwm=motif2matrix()

    simres= sim.scores(ncol(pwm),1000000)

    result=score.dist()

    df=data.frame(score=c(simres[[1]],result[[1]]),
                    prob=c(simres[[2]],result[[2]]),
                    model=c(rep("empiric",length(result[[1]])),
                    rep("dynprog",length(result[[1]]))))

    postscript(file=sprintf("%s/scoredist_account_%s_d%d.eps",figdir, pwmname,
                                                    variants$markov[var]),
                                            height=4,width=6)
    pl=ggplot(df,aes(x=score,y=prob,col=model))+
    labs(y=expression("P"["B"]*"(Score)"),x="Score")+
    geom_line()+
    theme_bw()+
            scale_color_manual(name="Distributions",
                    limits=c("empiric","dynprog"), 
                    values=c("red","black"),
                    labels=c(scorelabnames1[var],scorelabnames2[var]))+

    theme(text=element_text(size=20),
                            legend.position=c(.7,.7))
    print(pl)
    dev.off()

    r=which(cumsum(prob=result[[2]])>=.99)[1]
    smin=result[[1]][r]


    postscript(file=sprintf("%s/scoredist_account_%s_d%d_part.eps",figdir, pwmname,
                      variants$markov[var]),
                                            height=4,width=6)
    pl=ggplot(subset(df,score>=smin),aes(x=score,y=prob,col=model))+
    labs(y=expression("P"*"(Score)"),x="Score")+
    geom_line()+
    theme_bw()+
            scale_color_manual(name="Distributions",
                    limits=c("empiric","dynprog"), 
                    values=c("red","black"),
                    labels=c(scorelabnames1[var],scorelabnames2[var]))+

    theme(text=element_text(size=20),
                            legend.position=c(.7,.7))
    print(pl)
    dev.off()


    reldiff[var]=sum(-result[[2]][r]+simres[[2]][r])/sum(result[[2]][r])

    #sample counts
    sc=sim.counts(seqlen, maxhits, 10000)
    #compute scores
    ov=overlap.prob()
    cpc=comp.pois(seqlen, maxhits, 20, ov)

    dfc=data.frame(hits=c(1:length(sc[[1]]),1:length(cpc$dist))-1,
                            prob=c(sc[[1]],cpc$dist),
                            model=factor(c(rep("emp",length(sc[[1]])),
			rep("cp",length(cpc$dist)))))

    postscript(file=sprintf("%s/cntdist_account_%s_d%d.eps",figdir,pwmname,variants$markov[var]),
                                            height=4,width=6)
    pl=ggplot(dfc,aes(x=hits,y=prob,col=model))+
    geom_point()+
            labs(y=expression("P(Num. of hits)"),x="Num. of hits")+
            theme_bw()+
            scale_color_manual(name="Distributions",
                               limits=c("emp","cp"),
                               values=c("red","black"),
                    labels=c(countlabnames1[var],countlabnames2[var]))+

                theme(text=element_text(size=20),
                                    legend.position=c(.7,.7))

    print(pl)
    dev.off()

}
