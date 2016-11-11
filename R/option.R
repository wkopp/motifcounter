mdistOption=function(alpha, gran=0.1, ncores=1) {
  dummy=.C("mdist_option", as.numeric(alpha), as.numeric(gran),as.integer(ncores))
}
