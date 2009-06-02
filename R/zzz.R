.First.lib <- function(lib, pkg) {

library.dynam("bipartite", pkg, lib)

vers <- paste(sessionInfo()$otherPkg$bipartite$Version,".",sep="")

cat(paste("----------------------------------------------------------\nThis is bipartite",vers,"\nFor latest additions type: ?bipartite.\nFor citation please type: citation(\"bipartite\").\nHave a nice time plotting and analysing two-mode networks.\n----------------------------------------------------------\n\n"))
}

.Last.lib <- function() {

  library.dynam.unload("bipartite", paste(.libPaths(),"/bipartite",sep=""))

}