as.one.mode <- function(web, fill=0){
    # helper function
    #turns 2-mode matrix into 1-mode matrix
    # output can be used with the sna-package and its various indices
    #
    # after conversion, the object is called a "graph" in sna
    o <- matrix(fill, nrow=sum(dim(web)), ncol=sum(dim(web)))
    o[1:nrow(web), (nrow(web)+1):ncol(o)] <- web
    o[(nrow(web)+1):nrow(o), 1:nrow(web)] <- t(web)
    colnames(o) <- rownames(o) <- c(rownames(web), colnames(web))
    o
}