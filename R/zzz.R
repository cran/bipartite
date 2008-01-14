.First.lib <- function(lib, pkg) {

library.dynam("bipartite",pkg,lib)

}

.Last.lib <- function() {

library.dynam.unload("bipartite", paste(.libPaths(),"/bipartite",sep=""))
}