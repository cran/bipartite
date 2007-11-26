\name{second.extinct}

\alias{second.extinct}

\title{Secondary extinctions in bipartite networks}

\description{
Calculates the consequences of removing a species from a bipartite network. With plotting and slope-estimation functionality.
}

\usage{
second.extinct(web, participant = "high", method = "abun", nrep = 10,
  details = FALSE)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher
  trophic level species (columns) and lower trophic level species (rows). Usually
  this will be number of pollinators on each species of plants or number of parasitoids
  on each species of prey.}
  \item{participant}{\dQuote{\code{high}} (default) or \dQuote{\code{low}} or \dQuote{\code{both}},
  depending if you want to exterminate higher trophic level species, lower trophic
  level species, or a randomly chosen species of both levels; partial matching.}
  \item{method}{Random deletion of a species (\dQuote{\code{random}}) or according
  to its abundance, with least abundant going extinct first (\dQuote{\code{abundance}}; default). }
  \item{nrep}{Number of replicates of extermination sequence. Will not be considered
  for method abundance.}
  \item{details}{Logical; returns details, i.e.~for each replicate the sequence of
  secondary extinctions. If set to \code{FALSE} (default), replicated runs will be averaged. }
}
\details{
The procedure of this function is simple. For example imagine the web to represent
a pollination web, in which pollinators die one by one. Set all entries of a  column
to zero, see how may rows are now also all-zero (i.e. species that are now not pollinated
any more), and count these as secondary extinctions of the primary exterminated pollinator. 

Internally, each extermination is achieved by a call to \code{\link{extinction}},
followed by a call to \code{\link{empty}}, which counts the number of all-zero columns
and rows.

Although written for pollination webs (hence the nomenclature of pollinator and plant),
it can be similarly applied to other types of bipartite networks. It is called by
\code{\link{networklevel}}.

\option{details=FALSE} will not work with \option{participant="both"}, since
sequences cannot be averaged, due to different lengths.
}

\value{
Function returns an object of class \dQuote{bipartite} to ensure proper working of
function plot.bipartite.

If \option{details=FALSE}, the returned object contains a matrix with columns representing
the number of species going extinct from one step to the next, averaged across all runs.

If \option{details=FALSE}, the returned object will be a list of matrices containing
the number of species going extinct at each step.

The objects attribute \dQuote{exterminated} gives the name of the trophic level (pollinator
or plant).
}

\references{
Memmott, J., Waser, N. M. and Price, M. V. 2004 Tolerance of pollination networks
to species extinctions. \emph{Proceedings of the Royal Society b} \bold{271}, 2605--2611.
 }

\author{ Carsten F. Dormann }

\note{
Note: The length of an extinction sequence is obviously given by the number of species
in the selected trophic level. When setting \option{participant="both"}, lengths will
differ for each replicate run, since it is unpredictable in which sequence species go
extinct, and hence how many secondary extinctions will pre-empt further primary extinctions.
}

\seealso{ \code{\link{networklevel}} calls \code{second.extinct}; \code{\link{extinction}}
and \code{\link{empty}} are internal helper functions, and \code{\link{slope.bipartite}}
calculates extinction slopes from the output of \code{second.extinct}.}

\examples{
data(Safariland)
(ex <- second.extinct(Safariland, participant="low", method="r", nrep=50, details=TRUE))
(ex <- second.extinct(Safariland, participant="low", method="r", nrep=50, details=FALSE))
}

\keyword{package}

