readMotif=function(pwm, pseudocount=0.01) {
  if (is.matrix(pwm)) {
    if (any(pwm<=0)) {
        stop("All entries in the matrix must be greater than zero")
    }
    if (nrow(pwm)!=4) {
        stop("The number of rows must be 4, representing the number nucleotides.")
    }
    pwm=pwm+pseudocount
    pwm=pwm/apply(pwm,2,sum)
    dummy=.C("mdist_loadmotif", as.numeric(pwm), nrow(pwm), ncol(pwm))
  } else if (is.character(pwm)) {
    sdummy=.C("mdist_motiffromfile", as.character(pwm),
        as.numeric(pseudocount))
  } else {
      stop("pwm must be a filename pointing that contains a PFM or a PFM matrix")
  }
}

deleteMotif=function() {
  dummy=.C("mdist_deleteMotif")
}

motif2matrix=function() {
	m=.Call("mdist_fetchMotif");
	return (m)
}

motifLength=function() {
  x=integer(1);
  res=.C("mdist_motiflength",x)
  return (res[[1]])
}

