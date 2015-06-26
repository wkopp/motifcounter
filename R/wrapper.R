.onLoad=function(libpath,pkgname) {
  #packageStartupMessage("loading ", libpath, '::',pkgname)
}

.onUnload=function(pkgpath) {
  #packageStartupMessage("unloading ", pkgpath )
  mdist.unload()
}

mdist.unload=function() {
# to be removed when the package is detached
  dummy=.C("RdestroyBackground")
  dummy=.C("Rdestroymotif")
  
  #to be removed alpha, beta, delta
}

#.on.exit=function() {
#}
