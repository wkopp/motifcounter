
scoreDist=function() {
  scorerange=integer(1)
  scorerange=.C("mdist_scorerange",as.integer(scorerange))[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  .C("mdist_scoredist",as.numeric(scores),as.numeric(dist))
}

scoreDistBf=function() {
  scorerange=integer(1)
  scorerange=.C("mdist_scorerange",as.integer(scorerange))[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  .C("mdist_scoredist_bf",as.numeric(scores),as.numeric(dist))
}

