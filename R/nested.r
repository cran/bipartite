nested <- function(web, method="binmatnest2", ..., rescale=FALSE){
  # a wrapper function to call any of the currently implemented measures of nestedness
  
#  \item{method}{One or more of the following: discrepancy, discrepancy2, binmatnest, binmatnest2, NODF, NODF2, C.score, checker, wine, ALL}
# \item{rescale}{Should the methods all be transformed, so that higher values mean higher nesting? Only in combination with method "ALL". Defaults to FALSE.}
#  \item{...}{Arguments passed on the respective method.}

  if (! any(method %in% c("binmatnest", "discrepancy", "binmatnest2", "discrepancy2", "NODF", "NODF2", "wine", "C.score", "checker", "ALL"))) stop("Typo? Unknown method!")
  if ("ALL" %in% method) index <- c("binmatnest", "discrepancy", "binmatnest2", "discrepancy2", "NODF", "NODF2", "wine", "C.score", "checker") else index <- method

  out <- NULL
  if ("binmatnest2" %in% index) out <- c(out, "binmatnest2"=nestedtemp(web)$statistic)

  if ("binmatnest" %in% index) out <- c(out, "binmatnest"=nestedness(web, ...)$temperature)
  
  if ("discrepancy2" %in% index) out <- c(out, "discrepancy2"=nesteddisc(web)$statistic)

  if ("discrepancy" %in% index) out <- c(out, "discrepancy"=unname(discrepancy(web)))
     
  if ("C.score" %in% index) out <- c(out, "C.score"=C.score(web, ...)) #normalised!
  
  if ("checker" %in% index) out <- c(out, "checker"=nestedchecker(web)$C.score) 
  # identical to C.score(., FALSE)
 
  if ("NODF2" %in% index) out <- c(out, "NODF2"=unname(nestednodf(web, order=TRUE)$statistic[3])) # the nestedness of the row/column-sorted matrix (probably makes more sense)

  if ("NODF" %in% index) out <- c(out, "NODF"=unname(nestednodf(web, order=FALSE)$statistic[3])) # the "original", as I had implemented it, too (in NODF)

  if ("wine" %in% index) out <- c(out, "wine"=wine(web, ...)$wine)
  
  if (rescale & ! "ALL" %in% method) warning("You requested rescaling, but you won't get it (unless you use method='ALL')!") 
  
  if (rescale & "ALL" %in% method) out <- abs(c(100,100,0,0,0,0,0,0,0 ) - out)
  
  out
  
}
# example:
#data(Safariland)
#nested(Safariland, "ALL")
#nested(Safariland, "ALL", rescale=TRUE)
#
#nested(Safariland, c("C.score", "checker"), normalise=FALSE)