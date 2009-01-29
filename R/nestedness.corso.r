#the corso nestedness

nestedness.corso <- function(web, weighted=FALSE, reps=500){
    #this function implements the nestedness calculation according to Corso et al. (2008, arXiv: 0803.0007v1 [physics.bio-ph] 29 Feb 2008)
    #
    # there are two great news to nestedness freaks: 1. it is fast; 2. it is indifferent to transposing of the data!
    #
    # we can also use a weighted option, where a cell is weighted according to the number of observations it contains (as yet incomplete)
  
   d.weighted <- function(web){
    # helper function to pack web, calculate distances:
        web.packed <- web[order(rowSums(web>0), decreasing=TRUE), order(colSums(web>0), decreasing=TRUE), drop=FALSE]
        # note: there is no way to deal with ties implemented as yet!
        # also, for weighted webs, it is the number of LINKS that determine theorting, not the number of interactions! (is that intuitive? at least, that's what they write)
        L1 <- nrow(web); L2 <- ncol(web)
        coords <- which(web.packed>0, arr.ind=TRUE)
        coords.square <- coords
        coords.square[,1] <- (coords[,1]-1)/L1 + 1/(2*L1)
        coords.square[,2] <- (coords[,2]-1)/L2 + 1/(2*L2)

        # calculate dependencies:
        Pc <- web.packed/matrix(rowSums(web.packed), ncol=ncol(web), nrow=nrow(web), byrow=FALSE)
        Pr <- web.packed/matrix(colSums(web.packed), ncol=ncol(web), nrow=nrow(web), byrow=TRUE)
        # calculate weighted distances:
        d <- sum(coords.square[,1]*Pr[coords] + coords.square[,2]*Pc[coords]) #/ (sum(web>0))
        d
    }

    
    web <- as.matrix(web)
    web <- empty(web)
    if (weighted==FALSE) web <- (web>0) * 1 # turns quantitative matrix into a binary matrix
    colnames(web) <- rownames(web) <- NULL # just to get rid of species names
    
    
    #1. pack the matrix:
    web.packed <- web[order(rowSums(web>0), decreasing=TRUE), order(colSums(web>0), decreasing=TRUE), drop=FALSE]
    
    #2. project into unit square:
    # in fact, we now turn the matrix into a data.frame of coordinates, but only for the ones (since the distance to an empty cell is 0):
    L1 <- nrow(web); L2 <- ncol(web)
    coords <- which(web.packed>0, arr.ind=TRUE)
    coords.square <- coords
    coords.square[,1] <- (coords[,1]-1)/L1 + 1/(2*L1)
    coords.square[,2] <- (coords[,2]-1)/L2 + 1/(2*L2)
    
    N <- sum(web>0)
    
    #3. calculate "Manhattan" distances for coords.square (I'm not sure about this Manhattan business. In the paper, they use simply the sum of coordinates, but Manhattan is defined differently (see vegdist). I stick to the paper description, because it make sense in this context.):
    if (weighted==FALSE) {
      d <- sum(rowSums(coords.square)) # yes, I know that I don't need to sum over rows first; still, that is the way it is described in the paper; if we, at some point, want to weight the distance be the number of observations (ha!), then we need to weight it after summing the rows
      
        #4. d.rand = sum(d) = N #read the paper!
        d.rand <- sum(web)
    }

    #5. pack the matrix maximally and calculate d.min: the detail here is to fill the square matrix, not the original!
 #  also notice that the number of interactions has NO say in the sequence, only the distance to centre
    all.cells <- which(!is.na(web), arr.ind=TRUE)
    all.cells.square <- all.cells
    all.cells.square[,1] <- (all.cells[,1]-1)/L1 + 1/(2*L1)
    all.cells.square[,2] <- (all.cells[,2]-1)/L2 + 1/(2*L2)
    filling.seq <- order(rowSums(all.cells.square), decreasing=FALSE)
    
    max.packed <- matrix(0, nrow=nrow(web), ncol=ncol(web))
    max.packed[all.cells[filling.seq[1:N],]] <- 1
    
    coords.max <- which(max.packed==1, arr.ind=TRUE)
    coords.max.square <- coords.max
    coords.max.square[,1] <- (coords.max[,1]-1)/L1 + 1/(2*L1)
    coords.max.square[,2] <- (coords.max[,2]-1)/L2 + 1/(2*L2)

    if (weighted==FALSE) { d.min <- sum(rowSums(coords.max.square)) }
    if (weighted){
        d <- d.weighted(web)  # gives d for the packed web
        
        max.packed[all.cells[filling.seq[1:N],]] <- sort(web[web>0], decreasing=TRUE)
        d.min <- d.weighted(max.packed)
   
        # "heuristic" search for maximally un-nested web:
        random.d <- sapply(shuffle.web(web, reps), d.weighted)
        d.rand <- quantile(random.d, 0.95) # The 0.95% quantile is more robust to lower levels of randomisations than the max. This will yield a consistent overestimation of the "true" nestedness. Until a cool idea for finding the real maximum is available, this is better, in my view, than the max (thanks to Bernd for this idea).
        
    } 
    
    
    #6. range value:
    eta <- (d-d.min)/(d.rand-d.min) #done
    unname(round(eta, 2)) # that is all accuracy there is!
}

#nestedness.corso(Safariland)
#nestedness.corso(Safariland, weighted=TRUE)