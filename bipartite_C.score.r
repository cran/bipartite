#C-score (Stone & Roberts 1990, as in Gotelli & Rohde 2002)
#require(vegan)
#?designdist

C.score <- function(web, normalise=T, FUN=mean, ...){
    # calculates the C-score for all pollinator species; the C-score represents
    # the average number of checkerboard units for each unique specis pair.
    # (Stone & Roberts 1990; here taken from Gotelli & Rohde 2002)
    # these data are extremely skewed! A mean is, hm, not exactly what I would
    # suggest, hence the FUN-option allows for other summaries (try, e.g., hist).
    # option normalise ranges the index between 0 (no complementarity) and 1 (perfect
    # distinctness).
    # ... to be passed on to FUN (If a species occurs on all sites, then its distance
    # to any other will be 0, and maxD will also be 0. Since the division by 0 will
    # lead to NAs, use the ellipses to pass on na.rm to mean or similar functions. )
    # Carsten F. Dormann, Dec. 2007
    web <- web>0 # this whole concept works only on binary data!
    D <- designdist(t(web), method="(A-J)*(B-J)", terms="minimum")
    # The minimum value for Ds is 0, for the special case were all species use the
    # hosts exactly co-occurringly.
    # The maximum value for Ds in each comparison is AB, when they are exactly
    # complementary and hence J=0. However, if (A+B)>length of vector(L), then there
    # will be some co-occurrences and hence J>0=(A+B-L). The general maximum
    # then becomes (A-A-B+L)(B-A-B+L)=(L-B)(L-A). For (A+B)<L, maximum is AB.

    if (normalise){
      L <- ncol(t(web))
      maxD <- designdist(t(web), method="ifelse(L<(A+B),(L-A)*(L-B), (A*B))", terms="minimum")
      D <- D/maxD
    }
    FUN(D, ...)
}
# example:
m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,T)
C.score(m)

data(Safariland)
C.score(Safariland) # checkerboardness of pollinators
C.score(t(Safariland), FUN=hist) # checkerboardness of plants as histogram


V.ratio <- function(web){
  # invented by Schluter (1984)
  # be careful with the interpretation!!!
  web <- web>0
  
  N <- ncol(web) #number of species
  M <- nrow(web) #number of "samples"
  Tj <- colSums(web) #species "density" per sample j
  ni <- rowSums(web) #nr. of individuals per species i

  sigma.i <- ni/N*(1-ni/N)
  ST.2 <- 1/N*sum((Tj-mean(Tj))^2)
  V <- ST.2 / sum( (sigma.i)^2)
  V
}
#example:
V.ratio(Safariland)



nulls <- null.equip(100, t(Safariland))
null.C <- sapply(nulls, C.score) # takes a while!
fivenum(null.C) #[1] 0.03324430 0.07816549 0.10791768 0.13446594 0.29363426
hist(null.C, xlim=c(0, 1))
abline(v=C.score(t(Safariland))) # 0.77, and competition again!!!!


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# better null models:
#
# 1. keep row-sums constant, but allow column sums to freely vary:


web <- Safariland


null.equip <- function(N, sourceweb, prob=NULL){
    # implements a col-sum=constant null model
    # Ecologically, this is equivalent to treating hosts as equivalent and
    # allocating pollinators equiprobably to them. The number of pollinators does
    # NOT change, but their allocation to hosts does.
    #
    # In consequence, some columns will be empty!
    # prob allows a weighting by abundance of plants
    web <- (sourceweb>0)*1
    attr(web, "dimnames") <- NULL # get rid of names
    out <- as.list(1:N)
    for (i in 1:N) out[[i]] <- apply(web, 2, sample, prob=prob)
    out
}
null.equip(1, Safariland, rowSums(Safariland))

nulls <- null.equip(100, t(Safariland), prob=rowSums(t(Safariland)))
null.C <- sapply(nulls, C.score, normalise=T, na.rm=T)
fivenum(null.C) #[1] 0.2734649 0.3869649 0.4340690 0.4736084 0.5606600
hist(null.C, xlim=c(0, 1))
abline(v=C.score(t(Safariland), normalise=T)) # 0.56, and: no competition !


# 2. the swap-null model
# keeps row- and column-totals constant ("fixed-fixed")
null.swap <- function(N, web, swaps=500){
    web <- (web>0) *1
    attr(web, "dimnames") <- NULL
    m <- matrix(c(0,1,1,0), 2, 2, T)
    out <- as.list(1:N)
    for (j in 1:N){
      newweb <- web
      counter <- 0
      while (counter < swaps) {
          rowid <- sample(nrow(newweb), 2)
          colid <- sample(ncol(newweb), 2)
          subs <- newweb[rowid, colid]
          if (identical(subs*1, m)){
             newweb[rowid, colid] <- matrix(c(1,0,0,1), 2, 2, byrow=T)
             counter <- counter + 1
             # cat(counter, " ")
          }
      }
      out[[j]] <- newweb
    }
    out
}
null.swap(2, Safariland)

nulls <- null.swap(100, Safariland)
null.C <- sapply(nulls, C.score, na.rm=T) # takes a while!
fivenum(null.C) #[1] 0.007834758 0.080246914 0.116096866 0.152184236 0.219373219
hist(null.C, xlim=c(0,1))
abline(v=C.score(Safariland)) # 0.5776353: far higher than any of the r2dtable!
# this means, there are more checkerboard combinations than would be expected by
# chance alone -->> COMPETITION among pollinators for plants!


# 3. correct for abundance of plants:
nulls <- r2dtable(100, r=rowSums(Safariland), c=colSums(Safariland))
null.C <- sapply(nulls, C.score, na.rm=T) # takes a while!
fivenum(null.C) #[1] 0.007834758 0.080246914 0.116096866 0.152184236 0.219373219
hist(null.C, xlim=c(0,1))
abline(v=C.score(Safariland)) # 0.5776353: far higher than any of the r2dtable!
