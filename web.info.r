web.info <- function(web, name=substitute(web), author=NULL, licence=NULL, raw=FALSE){
  # generates information for the analysed network.
  # When given this information, null models can be compiled and compared to the
  # output of "networklevel".
  #
  # The idea is that for someone unwilling/unable to share the web-data itself can still share the analyses (networklevel, specieslevel, etc.) and these can be compared to null models.
  #
  # The information returned from this function will STILL allow anyone to reconstruct the network itself!! Marginal sums and absolute entries suffice! But: all names are purged. 
  # Web entries are sorted to prevent a reconstruction of the web.
  #
  # web     a bipartite network matrix
  # name    name of the network; defaults to the name used when calling the function
  # author  optional name of contact (possibly with email address)
  # licence optional statement of copyright licence (e.g. "GPL")
  # raw     also include the original, raw data; defaults to FALSE.
  
  rs <- unname(rowSums(web))
  cs <- unname(colSums(web))
  w <- as.vector(web)
  entries <- sort(w[w>0], decreasing=TRUE)
  out <- list("name"=name, "rowsums"=rs, "colsums"=cs, "entries"=entries)
  
  if (!is.null(author)) out["author"] <- author
  if (!is.null(licence)) out["licence"] <- licence
  if (raw) out[["raw data"]] <- web

  return(out)
}


data(Safariland)
web.info(Safariland, author="Diego Vazquez", licence="GPL", raw=TRUE)
