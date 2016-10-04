
r=3
paperdir="/project/Tf2motif/home/comppois_man/bioinfo/figures"

write(
paste(
markov, sprintf("$10^%d$",lalpha), seqlen[1],
#signif(sum(abs(pe$prob-pk$prob)),digits=r),
#signif(sum(abs(pe$prob-pp$prob)),digits=r),
#signif(sum(abs(pe$prob-pb$prob)),digits=r),
signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),sep=" & "),
file=sprintf("%s/%s_10-%d_m%d_5p.tabpart", 
						 paperdir,pwmname,lalpha,markov))

write(
paste(
markov, sprintf("$10^%d$",lalpha), seqlen[1],
signif(sum(abs(pe$prob-pk$prob)),digits=r),
signif(sum(abs(pe$prob-pp$prob)),digits=r),
signif(sum(abs(pe$prob-pb$prob)),digits=r),
#signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
#signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
#signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),
sep=" & "),
file=sprintf("%s/%s_10-%d_m%d_tv.tabpart", 
						 paperdir,pwmname,lalpha,markov))

average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

