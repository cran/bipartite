extinction <- function(web, participant="both", method="random", ext.row=NULL, ext.col=NULL){
  # function to simulate extinctions of pollination web participants
  # web contains a matrix (rows=plants, cols=pollinators)
  # participant selects which group is to suffer extinctions: "plant", "pollinator" or "both";
  #        partial matching of strings allowed
  # method is either "random" (randomly a plant species is exterminated) or
  #        "abundance" (plant species with the lowest number of interactions
  #        is exterminated first; idea is that interaction number is proportional
  #        to true abundance); 
  #		   "degree", eliminating the most highly connected species first; or
  #		   "external", in which case ext.row and ext.col need to give the sequence by which species are to be eliminated (e.g. by inverse of independently measured abundance) ; partial matching of strings allowed
  # based on how I remember the idea of Jane Memmott's paper; C.F. Dormann, 6 Mar 2007

  # Error checking for partial matching of options:
  partis <- c("lower", "higher", "both")
  partis.match <- pmatch(participant, partis)
  if (is.na(partis.match)) stop("Choose participant: lower/higher/both.\n")

  meths <- c("random", "abundance", "degree", "external")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) stop("Choose extinction method: random/abundance/degree.\n")

  nr <- NROW(web); nc <- NCOL(web)

  if (partis.match==3) partis.match <- sample(2,1) #randomly draw a participant

  # construct row & column indices for random/abundance removal
  if (meths.match == 1) # random removal:
  {
      rexcl <- sample(nr,1) #removes a plant
      cexcl <- sample(nc,1) #removes a pollinator
      if (partis.match==1) web[rexcl, ] <- 0
      if (partis.match==2) web[, cexcl] <- 0

  } 
  if (meths.match == 2)               # removal by abundance:
  {
      rseq <- order(rowSums(web))
      cseq <- order(colSums(web))
      if (partis.match==1) web[rseq[1], ] <- 0
      if (partis.match==2) web[, cseq[1]] <- 0

  }
  if (meths.match == 3){     # removal by degree
      if (partis.match == 1){
          sequ <- rowSums(web>0)   
          which.ex <- which(sequ == max(sequ))
          # if two or more species have the same degree, draw one randomly:
          if (length(which.ex) > 1){ ex <- sample(which.ex, size=1)}  else {ex <- which.ex}
          web[ex, ] <- 0
      }
      if (partis.match == 2){
          sequ <- colSums(web>0)
          which.ex <- which(sequ == max(sequ))
          if (length(which.ex) > 1) ex <- sample(which.ex, size=1)  else ex <- which.ex
          web[, ex] <- 0
      }
     
  }
  if (meths.match == 4)               # removal by external sequence vector
  {
      rseq <- ext.row
      cseq <- ext.col
      if (partis.match==1) web[rseq[1], ] <- 0
      if (partis.match==2) web[, cseq[1]] <- 0

  }
  

  return(web)
}


#extinction(Safariland, "higher", "abun")
#extinction(Safariland, "higher", "deg")

