\name{degreedistr}
\alias{degreedistr}

\title{Fits functions to cumulative degree distributions of both trophic levels of a network.}

\description{
This function first calculates degrees for each species, then constructs a cumulative distribution with them, and finally fits three different functions to these distributions: exponential, power law and truncated power law. Coefficients and fits are returned.
}

\usage{
degreedistr(web, plot.it=TRUE, pure.call=TRUE, silent=TRUE)
}

\arguments{
  \item{web}{A bipartite network matrix.}
  \item{plot.it}{Logical; returns graphs of fits when set to TRUE (default).}
  \item{pure.call}{Logical; adjusts par for two panels (for TRUE) or leaves this to the wrapper function (FALSE).}
  \item{silent}{Logical; suppresses error reporting in the \code{try}-function around \code{nls}; defaults to TRUE.}
}

\details{
Jordano et al. (2003) proposed that plant-animal networks may show scale invariance, as indicated by the presence of a power law in species degrees. They report on consistently better fits of the truncated power law, hypothesising that such patterns may arise from morphological mismatch or phenological uncoupling.
Most problematic with the use of this particular approach is the extreme demand for data. The example web \code{Safariland} in this package is large (1130 interactions), but it provides only 5 different degree levels (for plants, only 4 for pollinators). Hence fitting three different non-linear functions to these few points is stretching it a bit.
}

\value{
 For both trophic levels, a table:
  \item{... trophic level dd fits}{Contains coefficient estimates, estimate's standard
      error and P-value, R2 and AIC for each of the three model fits, for the respective trophic level.}
}

\references{
Jordano, P., Bascompte, J. and Olesen, J. M. (2003) Invariant properties in coevolutionary networks of plant-animal interactions. \emph{Ecology Letters} \bold{6}, 69-81  }

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de} }

\note{
The truncated power law fits two coefficients: slope and cut-off. The function only returns the slope.
R2-values for non-linear fits are not well liked among statisticians! See the discussion the R-help list (e.g. \url{http://finzi.psych.upenn.edu/R/Rhelp02a/archive/90168.html}).
Finally, often data are too few to yield any fit. In this case the error message \dQuote{singular gradient} is returned to signalise this problem!

Post finally, yes, I am aware that degrees are integers and unlikely to be normally distributed, and that thus the \code{nls} procedure is not really a good idea. My (poor) excuse: I followed the implementation of the above-cited paper and do not believe enough in degree distributions (and power laws, for that matter) to implement a proper likelihood-based approach. Check out the \pkg{statnet} bundle for alternative approaches to this problem.
}

\seealso{ \code{\link{networklevel}}, where \code{degreedistr} is called (without picturing the results) }

\examples{
data(Safariland)
degreedistr(Safariland)
}

\keyword{ htest }

