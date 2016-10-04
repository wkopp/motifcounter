paperdir="project/Tf2motif/home/thesis/figures/tab_chapter2"

write(
paste(
markov, sprintf("$10^%d$",lalpha),
signif(sum(abs(pe$prob-pk$prob)),digits=r),
signif(sum(abs(pe$prob-pp$prob)),digits=r),
signif(sum(abs(pe$prob-pb$prob)),digits=r),
signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),sep=" & "),
file=sprintf("%s/%s_10-%d_m%d.tabpart", 
						 paperdir,pwmname,lalpha,markov))

average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

