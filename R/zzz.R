.First.lib <- function(lib, pkg) {

library.dynam("bipartite",pkg,lib)
cat("This is bipartite. Have a nice time plotting and analysing two-mode networks.")
}

.Last.lib <- function() {

library.dynam.unload("bipartite", paste(.libPaths(),"/bipartite",sep=""))

}