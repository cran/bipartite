`shuffle.web` <-
function(web, N){
  # shuffles entries of an interaction web, thereby maintaining connectence, but
  # changing marginal sums
  # to maintain dimensionality, interactions are first allocated to the diagonal
  # (thereby having all species in the web), then to the rest
  shuffle <- function(web){
      web <- as.matrix(web)
      web <- empty(web)
      dimdiff <- dim(web)[1]-dim(web)[2]
      # more columns than rows, i.e. dimdiff is negative
      if ((length(diag(web))+dimdiff) > sum(as.vector(web)>0))  stop("Too few entries in the web: less interactions than length of web diagonal.")
      out <- web
      out[,] <- 0
      shuf <- sample(as.vector(web)) #shuffle web entries
      nozero.index <- which(shuf!=0) #pick non-zeros
      diag(out) <- shuf[nozero.index[1:length(diag(out))]] #allocate to diagonal of new web
    
      # find remaining cols/rows beyond the diagonal:
      if (dimdiff < 0) colindex <- ((length(diag(out))+1) : dim(web)[2]) else colindex <- 1:dim(web)[2]
      if (dimdiff > 0) rowindex <- ((length(diag(out))+1) : dim(web)[1]) else rowindex <- 1:dim(web)[1]
    
      if (dimdiff < 0) rowposition <- sample(rowindex, length(colindex), replace=TRUE)
      if (dimdiff > 0) colposition <- sample(colindex, length(rowindex), replace=TRUE)
    
      # slot in non-zeros:
      if (dimdiff < 0) {
          for (i in 1:abs(dimdiff)){
            out[rowposition[i], colindex[i]] <- shuf[nozero.index[(length(diag(web))+i)]]
          }
      }
      if (dimdiff > 0){
          for (j in 1:abs(dimdiff)){
           out[rowindex[j], colposition[j]] <- shuf[nozero.index[(length(diag(web))+j)]]
          }
      }
    
      # allocate remaining entries:
      remains <- shuf[!(shuf %in% out)]  
#      remains <- shuf[nozero.index[(length(diag(web))+abs(diff(dim(out)))+1):length(nozero.index)]]
      option <- which(out == 0, arr.ind=TRUE)
      out[option[sample(dim(option)[1], length(remains)),]] <- remains
      colnames(out) <- rownames(out) <- NULL
      out
  }
  replicate(N, shuffle(web), simplify=FALSE)
}

