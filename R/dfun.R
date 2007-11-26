`dfun` <-
function(web, abundances=NULL){
    # The abundance vector allows to incorporate independent estimates of the
    # abundances of the HIGHER trophic level. In a pollination web, pollinator abundances
    # may be very different from those estimated by the interaction matrix column sums.
    # This has also, obviously, large consequences for the specialisation: A plant
    # being pollinated by a bee that is common on this plant, but very rare in general,
    # will show a low specialisation unless bee abundances are corrected for.
    # Data given in the abundance vector are here used in replacement for the row sums,
    # both in the d-function itself, as well as in the calculation of the minimum ds.

    web <- empty(web)
    #-- -- -- -- -- -- --
    # d uncorrected:
    d <- function(web, abundances=NULL){#
      #calculates the uncorrected d for the LOWER trophic level of a web

      Pprime <- web/rowSums(web, na.rm=TRUE) #observed proportions of each higher species on each lower species
      
      if (is.null(abundances)){
          #expected proportion of higher trophic level attributed to each species (availability of higher trophic level):
          Q <- colSums(web, na.rm=TRUE)/sum(web)
      } else {
          # if independent abundances for lower trophic level are given:
          Q <- abundances/sum(abundances)
      }

      Qmat <- matrix(rep(Q, NROW(web)), nrow=NROW(web), byrow=TRUE)
      rowSums(Pprime * log(Pprime/Qmat), na.rm=TRUE)  #Spezialisierung von Art i; formula is KL-distance

    }
    duncorr <- d(web, abundances=abundances)
    #-- -- -- -- -- -- --
    # d min: (Variation on Jochen's way)
      if (is.null(abundances)){
          #same number of interactions as web, but now as expected from marginal sums:
          exexpec <- outer(rowSums(web), colSums(web)/sum(web))
      } else {exexpec <- outer(rowSums(web), abundances/sum(abundances))}
      expec <- floor(exexpec)           # down-rounded expectation values

# Jochen's version (from H2max)
#      rs <- rowSums(web)
#      cs <- colSums(web)
#      newweb <- expec  # start new web
#      webfull <- matrix("nö", nrow(web), ncol(web)) # makes boolean web, set to 0
#      while (sum(newweb)<sum(web)) {
#         difexp <- exexpec-expec               #use for allocation ranking
#         webfull[which(rowSums(newweb)==rs),]="yo" # sets columns/rows with correct cs/rs to 1
#         webfull[,which(colSums(newweb)==cs)]="yo"
#         OK <- webfull=="nö" # matrix of potential cells
#         smallestpossible <- newweb==min(newweb[OK])  # find cell with lowest number of interactions (e.g. 0)
#         greatestdif <- max(difexp[smallestpossible & OK]) # find cell value with largest different between "is" and "should"
#         bestone <- which(OK & smallestpossible & difexp==greatestdif ) # find cell for all three conditions
#         if (length(bestone)>1) bestone <- sample(bestone,1) # select randomly a cell, if different are possible
#         newweb[bestone] <- newweb[bestone]+1 # put an interaction into that cell
#      }
#
      restuse <- sum(web) - sum(expec)  #number of interactions left to allocate
      #find empty cols/rows:
      emptycols <- which(colSums(expec)==0)
      emptyrows <- which(rowSums(expec)==0)
      if ((length(emptycols)+length(emptyrows))>restuse) warning("There are more empty cols/rows than spare interactions.")
      difexp <- exexpec-expec               #use for allocation ranking
      # put 1s into empty cols/rows (to keep the same number of species as in web):
      replaceindexC <- apply(as.matrix(difexp[,emptycols]), 2, which.max)
      expec[replaceindexC, emptycols] <- 1
      replaceindexR <- apply(as.matrix(difexp[emptyrows,]), 1, which.max)
      expec[emptyrows, replaceindexR] <- 1
      restuse <- sum(web) - sum(expec)  #number of interactions left to allocate
      difexp <- exexpec-expec               #use for allocation ranking
#      if (restuse>0) {
#          # this is not very elegant: it now allocates more interactions than it
#          # actually has (but at least returns a warning)
#          replaceindex <- which(difexp %in% sort(difexp, decreasing=TRUE)[1:restuse])
#          expec[replaceindex[1:restuse]] <- expec[replaceindex[1:restuse]] + 1  # this is the Dmin matrix
#      }
#
      rs <- rowSums(web); cs <- colSums(web)
      newweb=expec
      webfull <- matrix("no", nrow(web), ncol(web)) # makes boolean web, set to 0
      while (restuse>0){
          replaceindex <- match(difexp, sort(difexp, decreasing=TRUE))[1]
          webfull[which(rowSums(newweb)==rs),]="yes" # sets columns/rows with correct cs/rs to 1
          webfull[,which(colSums(newweb)==cs)]="yes"
          OK <- webfull=="no" # matrix of potential cells
          smallestpossible <- (expec==min(expec[OK]))  # find cell with lowest number of interactions (e.g. 0)
          greatestdif <- max(difexp[smallestpossible & OK]) # find cell value with largest different between "is" and "should"
          bestone <- which(OK & smallestpossible & difexp==greatestdif ) # find cell for all three conditions
          if (length(bestone)>1) bestone <- sample(bestone,1) # select randomly a cell, if different are possible
          newweb[bestone] <- newweb[bestone]+1 # put an interaction into that cell
          restuse <- restuse - 1
      }
      
      #expec <- expec[rowSums(expec)!=0, colSums(expec)!=0]
      #newweb=expec
      dmin <- ifelse(duncorr<d(newweb), duncorr, d(newweb))
     # if (any(duncorr < d(newweb))) warning("This is one of the extremely rare occassions when duncorr is for rounding reasons smaller than dmin (not to worry).")

    #-- -- -- -- -- -- --
    # d max:
    dmax <- log(sum(web)/rowSums(web))
    #ifelse(duncorr>log(sum(web)/rowSums(web)), duncorr, log(sum(web)/rowSums(web)))

    # This is NOT solved satisfactorily!!
    # Problem here is that dmax for cases when abundances are given are calculated
    # very differently to simply log(sum(web)/rowSums(web)) [as given in Blüthgen et al. 2006].
    # I simply do not know how to get to the correct maximum d-values for the different species!

    #-- -- -- -- -- -- --
    # d max:
    dcorr <- (duncorr-dmin)/(dmax-dmin)
    list("dprime"=dcorr, "d"=duncorr, "dmin"=dmin, "dmax"=dmax)
}

# dfun(Safariland, abundances=runif(ncol(Safariland)))
