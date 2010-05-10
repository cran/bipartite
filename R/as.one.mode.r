as.one.mode <- function(web, fill=0, project="full", weighted=TRUE){
    # helper function
    #turns 2-mode matrix into 1-mode matrix
    # output can be used with the sna-package and its various indices
    #
    # after conversion, the object is called a "graph" in sna

    projection <- function(web, weighted=TRUE){
      # projection helper function
      N <- NCOL(web)
      as.one.mode.web <- matrix(0, N, N)
      colnames(as.one.mode.web) <- rownames(as.one.mode.web) <- colnames(web)
      for (i in 1:N){
        for (j in 1:N){
          if (j >= i) next;
          set <- web[,c(i,j)]
          Ints <- pmin(set[,1], set[,2])
          Links <- sum(Ints > 0)
          Weight <- sum(Ints)
          if (weighted) W <- Weight else W <- 1
          if (Links > 0) as.one.mode.web[i,j] <- W
        }
      }
      as.one.mode.web <- as.one.mode.web + t(as.one.mode.web)
      as.one.mode.web
    }


    # maintain full information, pretend all species could interact:
    if (project == "full"){
      if (!weighted) web <- web>0
      o <- matrix(fill, nrow=sum(dim(web)), ncol=sum(dim(web)))
      o[1:nrow(web), (nrow(web)+1):ncol(o)] <- as.matrix(web)
      o[(nrow(web)+1):nrow(o), 1:nrow(web)] <- t(web)
      colnames(o) <- rownames(o) <- c(rownames(web), colnames(web))
    }
    
    # project to one mode for higher trophic level
    if (project == "higher"){
      o <- projection(web, weighted=weighted)
    }
    
    # project to one model for lower trophic level
    if (project == "lower"){                                                
      o <- projection(t(web), weighted=weighted)
    }   
    o
}


#m <- matrix(0,9,7) # Dalgaard's example
#m[1, 1:4] <- 1
#m[2, c(3, 5:7)] <- 1
#m[3, 7] <- m[4,1] <- m[5,2] <- m[6,4] <- m[7,3:4] <- m[8,6] <- m[9, 4:5] <- 1
#m                             
#
#m.1 <- projection(m)

