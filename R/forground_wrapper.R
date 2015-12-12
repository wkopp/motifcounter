motif.length=function() {
  x=integer(1);
  res=.C("Rmotiflength",x)
  return (res[[1]])
}

read.motif=function(pwmfile, format, pseudocount=0.01) {
  sdummy=.C("Rmotiffromfile", as.character(pwmfile),
    as.numeric(pseudocount), as.character(format))
}

matrix2motif=function(data) {
	if (!is.matrix(data)) {
		stop("data must be a matrix")
	}
	data=data/apply(data,2,sum)
  dummy=.C("Rloadmotif", as.numeric(data), nrow(data), ncol(data))
  #if (is.matrix(data) && all(apply(data,2,sum)==rep(1,ncol(data)))) {
  #} else {
  #  warning("provided data must be a 4xM matrix, with the columns summing to 1.")
  #}
}

delete.motif=function() {
  dummy=.C("Rdestroymotif")
}

motif2matrix=function() {
	m=.Call("fetchMotif");
	return (m)
}
