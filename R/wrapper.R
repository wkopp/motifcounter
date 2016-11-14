.onLoad=function(libpath,pkgname) {
  #packageStartupMessage("loading ", libpath, '::',pkgname)
}

.onUnload=function(pkgpath) {
  #packageStartupMessage("unloading ", pkgpath )
#  mdist.unload()
}

