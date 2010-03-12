\encoding{latin1}
\name{ND}

\alias{ND}
\alias{BC}
\alias{CC}

\title{Normalised degree, betweenness and closeness centrality}

\description{
 Calculates normalised degrees, and two measures of centrality, betweenness and closeness. These two are based on one-mode representations of the network and invoke functions from \pkg{sna}.
}

\usage{
ND(web, normalised=TRUE)
BC(web, rescale=TRUE, ...)
CC(web, cmode="suminvdir", rescale=TRUE, ...)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species
    as columns and number of interactions as entries.}
  \item{normalised}{Shall the degrees be normalised? If so (default), the degree for a species is divided by the number of potential partners-1 (see, e.g., Mart�n Gonz�lez et al. 2009).}
  \item{rescale}{If TRUE (default), centrality scores are rescaled such that they sum to 1.}
  \item{cmode}{String indicating the type of betweenness centrality being computed (directed or undirected geodesics, or a variant form - see help for \code{closeness} in \pkg{sna} for details). The default, \option{"suminvdir"}, uses a formula that can also be applied to disconnected (=compartmented) graphs. Other cmodes cannot.}
  \item{...}{Options passed on to \code{betweenness} and \code{closeness}, respectively.}
}

\details{
  These functions are convinience functions to enable easy reproduction of the type of analyses by Mart�n Gonz�lez et al. (2009). BC and CC are wrappers calling two functions from \pkg{sna}, which uses one-mode, rather than bipartite data. 
}

\value{
  A list with two entries, ``lower'' and ``higher'', which contain a named vector of normalised degrees, betweenness centrality and closeness centrality, respectively. The lower-entry contains the lower trophic level species, the higher analogously the higher trophic level species.
}

\references{
 Mart�n Gonz�les, A.M., Dalsgaard, B. and Olesen, J.M. 2009. Centrality measures and the importance of generalist species in pollination networks. \emph{Ecological Complexity}, in press (doi:10.1016/j.ecocom.2009.03.008)
}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de} }

\note{ 
Experimental. Should work most of the time, but not necessarily always. Also, on trials with the same data as those of Mart�n Gonz�lez et al. (2009), numerical values differed slightly. Whether this is due to rounding errors, different non-linear least square fits in JMP and R or whatever I cannot tell. See example for my attempt to reproduce their values for the network ``Azores'' (aka \code{\link{olesen2002flores}}).
}

\seealso{
  \code{centralization}, \code{betweenness} and \code{closeness} in \pkg{sna}; \code{\link{specieslevel}} which calls them
}

\examples{
## example:
data(olesen2002flores)
(ndi <- ND(olesen2002flores))
(cci <- CC(olesen2002flores))
(bci <- BC(olesen2002flores))

cor.test(bci[[1]], ndi[[1]], method="spear") # 0.779
cor.test(cci[[1]], ndi[[1]], method="spear") # 0.826

cor.test(bci[[2]], ndi[[2]], method="spear") # 0.992
cor.test(cci[[2]], ndi[[2]], method="spear") # 0.919

## PLANTS:
bc <- bci[[1]]
cc <- cci[[1]]
nd <- ndi[[1]]
# CC:
summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) # lower RSE
summary(nls(cc ~ c*nd^d, start=list(c=0.072,d=0.2))) 
# BC:
summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE

## ANIMALS:
bc <- bci[[2]]
cc <- cci[[2]]
nd <- ndi[[2]]
# CC:
summary(nls(cc ~ a*nd+b, start=list(a=1,b=1)))  # lower RSE
summary(nls(cc ~ c*nd^d, start=list(c=0.2,d=2))) 
# BC:
summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
}

\keyword{package}

