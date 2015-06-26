
#interface functions for motif representation
read.background=function(file, order=1) {
  if (regexpr(pattern=".fasta$",file)>=0) {
    dummy=.C("Rmakebg", as.character(file), as.integer(order))
  } else {
    dummy=.C("RreloadBackground", as.character(file))
  }
}

store.background=function(file) {
  dummy=.C("Rstorebg",as.character(file));
}

print.background=function() {
  dummy=.C("RprintBackground");
}

delete.background=function() {
  dummy=.C("RdestroyBackground")
}
