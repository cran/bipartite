\encoding{latin1}

\name{shuffle.web}
\alias{shuffle.web}

\title{ Shuffle web entries  }
\description{
Shuffles (= relocates) entries in a web matrix whilst maintaining the dimensionality.
}
\usage{
shuffle.web(web, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{web}{ An interaction matrix.}
  \item{N}{ Number of desired shuffled matrices. }
}
\details{
 This function is designed to behave similar to \code{r2dtable}, i.e. it returns a list of randomised matrices. 
 In contrast to \code{r2dtable} is does not keep marginal sums constant! 
 
 This function is thought of as a nullmodel for the analysis of bipartite webs. It keeps two web properties
 constant: The number of interactions and the number of links (and hence connectance). A comparison of
 \code{shuffle.web}- and \code{r2dtable}-based webs allows to elucidate the effect of marginal sums.
}
\value{
 A list of N randomised matrices with the same dimensions as the initial web.
}

\references{
 There is no reference for this specific nullmodel. For a very nice and thorough overview of nullmodels in general see

 Gotelli, N. J., and G. R. Graves. 1996. Null Models in Ecology. Smithsonian Institution Press, Washington D.C.
 
 For nullmodels and their application to webs/networks see, e.g.:

 Vázquez, D. P., and M. A. Aizen. 2003. Null model analyses of specialization in plant-pollinator interactions. Ecology 84:2493-2501.

 Vázquez, D. P., C. J. Melián, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecological networks. Oikos 116:1120-1127.

}

\author{ Carsten F. Dormann <carsten.dormann@ufz.de> }

\note{ \code{shuffle.web} is not an extremely intelligent nullmodel. You may want to think of a better one for your specific application!
}

\seealso{ \code{\link{r2dtable}} }

\examples{

data(Safariland)

shuffle.web(Safariland, N=2)

}
\keyword{ package }

