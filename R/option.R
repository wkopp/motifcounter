mdist.option=function(alpha, gran,ncores=1) {
  dummy=.C("Roption", as.numeric(alpha), as.numeric(gran),as.integer(ncores))
}
