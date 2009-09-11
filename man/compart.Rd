\name{compart}
\alias{compart}

\title{ Detects compartments }
\description{
 Finds number of compartments, based on multivariate ordination techniques, and labels interactions according to the compartment they belong to.
}

\usage{
compart(web)
}

\arguments{
  \item{web}{A bipartite interaction web, i.e.~a matrix with higher (cols) and lower (rows) trophic levels.}
}

\details{
Internal function, to be called by \code{\link{networklevel}}. 
}

\note{
Note that up to (and including) version 0.85 we used a code based on correspondence analysis (see Lewinsohn et al. 2006). This is, however, faulty for webs with many same-linked species. Hence we resort to a brute-force search for compartments, which is orders of magnitude slower, but at least works correctly.

For large networks, the recursive call of compartment my cause an error. Use \code{options(expressions=50000)} at the beginning of the script to allow for deeper nesting.
}

\value{
Returns a list with two entries:
  \item{cweb}{A matrix similar to \code{web}, but now with compartment numbers instead of interaction values.}
  \item{ncompart}{The number of compartments.}
}

\references{
Lewinsohn, T. M., P. I. Prado, P. Jordano, J. Bascompte, and J. M. Olesen (2006) Structure in plant-animal interaction assemblages. \emph{Oikos} \bold{113}, 174--184
}

\author{ Carsten F. Dormann <carsten.dormann@ufz.de> }

\seealso{ See also \code{\link{networklevel}}. }

\examples{
# make a nicely compartet web:
web <- matrix(0, 10,10)
web[1,1:3] <- 1 
web[2,4:5] <- 1 
web[3:7, 6:8] <- 1
web[8:10, 9:10] <- 1
web <- web[-c(4:5),] #oh, and make it asymmetric!
web <- web[,c(1:5, 9,10, 6:8)] #oh, and make it non-diagonal
compart(web)

# or, standard, use Safariland as example:
data(Safariland)
compart(Safariland)
}
\keyword{ package}

