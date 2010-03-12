# normalised degree
ND <- function(web, normalised=TRUE){
    # calculates the (normalised) degree of a species
    # by Carsten F. Dormann, 14 Dec 2010
    web <- (web > 0) * 1
    k <- sum(web)
    dlower <- rowSums(web)
    dhigher <- colSums(web)  
    Nlow <- Nhigh <- 2 # effectively unnormalised
    if (normalised){
      Nlow <- length(dhigher)+1
      Nhigh <- length(dlower)+1
    }
    low <- dlower/(Nlow-1); names(low) <- rownames(web)
    high <- dhigher/(Nhigh-1); names(high) <- colnames(web)
    list("lower"=low, "higher"=high)
}

CC <- function(web, cmode="suminvdir", rescale=TRUE, ...){
    # closeness centrality as used in Gonzales et al. (2009)
    require(sna)
    # by Carsten F. Dormann, 14 Dec 2010
    # ... options passed on to closeness in package sna
    # uses a version to calculate centrality that allows for disconnected graphs
    w <- as.one.mode(web)
    cc <- sna::closeness(w, cmode=cmode, rescale=rescale, ...)
    low <- cc[1:nrow(web)]; names(low) <- rownames(web)
    high <- cc[(nrow(web)+1):length(cc)]; names(high) <- colnames(web)
    list("lower"=low, "higher"=high)
}

BC <- function(web, rescale=TRUE, ...){
    # betweenness centrality as used in Gonzales et al. (2009)
    require(sna)
    # by Carsten F. Dormann, 14 Dec 2010
    # ... options passed on to closeness in package sna; particularly: rescale=TRUE!
    w <- as.one.mode(web)
    bc <- sna::betweenness(w, rescale=rescale, ...)
    low <- bc[1:nrow(web)]; names(low) <- rownames(web)
    high <- bc[(nrow(web)+1):length(bc)]; names(high) <- colnames(web)
    list("lower"=low, "higher"=high)

}



## example:
#require(bipartite)
#data(olesen2002flores)
#(ndi <- ND(olesen2002flores))
#(cci <- CC(olesen2002flores))
#(bci <- BC(olesen2002flores))
#
#cor.test(bci[[1]], ndi[[1]], method="spear") # 0.779
#cor.test(cci[[1]], ndi[[1]], method="spear") # 0.826
#
#cor.test(bci[[2]], ndi[[2]], method="spear") # 0.992
#cor.test(cci[[2]], ndi[[2]], method="spear") # 0.919
#
### PLANTS:
#bc <- bci[[1]]
#cc <- cci[[1]]
#nd <- ndi[[1]]
## CC:
#summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) # lower RSE
#summary(nls(cc ~ c*nd^d, start=list(c=0.02,d=2))) 
## BC:
#summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
#summary(nls(bc ~ c*nd^d, start=list(c=0.02,d=2))) # lower RSE
#
### ANIMALS:
#bc <- bci[[2]]
#cc <- cci[[2]]
#nd <- ndi[[2]]
## CC:
#summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) 
#summary(nls(cc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
## BC:
#summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
#summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
#
#
##see also, for whole web measures:
#centralization(web, "degree")
#centralization(web, "betweenness")