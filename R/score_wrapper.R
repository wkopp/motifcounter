
scoreDist=function() {
  scorerange=integer(1)
  scorerange=.C("Rscorerange",as.integer(scorerange))[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  .C("Rscoredist",as.numeric(scores),as.numeric(dist))
}

