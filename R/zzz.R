.First.lib <- function(lib, pkg) {

library.dynam("bipartite",pkg,lib)
cat("----------------------------------------------------------\nThis is bipartite.\nFor citation please type citation(\"bipartite\").\nHave a nice time plotting and analysing two-mode networks.\n----------------------------------------------------------\n\n")
}

.Last.lib <- function() {

library.dynam.unload("bipartite", paste(.libPaths(),"/bipartite",sep=""))

}