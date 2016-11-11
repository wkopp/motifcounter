.onLoad=function(libpath,pkgname) {
  #packageStartupMessage("loading ", libpath, '::',pkgname)
}

.onUnload=function(pkgpath) {
  #packageStartupMessage("unloading ", pkgpath )
  mdist.unload()
}

mdist.unload=function() {
# to be removed when the package is detached
  dummy=.C("mdist_deleteBackground")
  dummy=.C("mdist_deleteMotif")
  
  #to be removed alpha, beta, delta
}

#.on.exit=function() {
#}
