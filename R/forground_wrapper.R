readMotif=function(pwm, pseudocount=0.01) {
  if (is.matrix(pwm)) {
    if (any(pwm<=0)) {
        stop("All entries in the matrix must be greater than zero")
    }
    if (nrow(pwm)!=4) {
        stop("The number of rows must be 4, representing the number nucleotides.")
    }
    pwm=pwm/apply(pwm,2,sum)
    dummy=.C("Rloadmotif", as.numeric(pwm), nrow(pwm), ncol(pwm))
  } else {
    sdummy=.C("Rmotiffromfile", as.character(pwm),
        as.numeric(pseudocount))
  }
}

matrix2motif=function(data) {
  if (!is.matrix(data)) {
     stop("data must be a matrix")
  }
}

deleteMotif=function() {
  dummy=.C("Rdestroymotif")
}

motif2matrix=function() {
	m=.Call("fetchMotif");
	return (m)
}

motifLength=function() {
  x=integer(1);
  res=.C("Rmotiflength",x)
  return (res[[1]])
}

