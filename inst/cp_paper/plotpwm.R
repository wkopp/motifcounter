postscript(file=sprintf("%s/%s_10-%d_m%d.eps", figdir,pwmname,lalpha,markov),
					 width=6,height=4)
pl=ggplot(df,aes(x=hits,y=prob, col=model, shape=model))+
  labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
  scale_colour_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c("black", "blue", "red", "grey"),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  scale_shape_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c(20,19,17,3),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
  theme(text=element_text(size=20),
  			legend.position=c(0.7,.7),
  			legend.text=element_text(size=18),
  			legend.key.size=unit(.4,"inches"))
print(pl)
dev.off()

postscript(file=sprintf("%s/%s_10-%d_m%d_part.eps", figdir,pwmname,lalpha,markov),
					 width=6,height=4)

plp=ggplot(dfp,aes(x=hits,y=prob, col=model, shape=model))+
  labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
  scale_colour_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c("black", "blue", "red", "grey"),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  scale_shape_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c(20,19,17,3),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
  theme(text=element_text(size=20),
  			legend.position=c(0.7,.7),
  			legend.text=element_text(size=18),
  			legend.key.size=unit(.4,"inches"))

print(plp)
dev.off()
r=3

