\encoding{latin1}
\name{dfun}

\alias{dfun}

\title{Calculates standardised specialisation index d' (d prime) for each species in the lower trophic level of a bipartite network.}
\description{
This function returns the specialisation index d', which expresses how specialised a given species is in relation to what resources are on offer.
}

\usage{
dfun(web, abundances=NULL)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{abundances}{A vector of abundances for the higher trophic level, usually from additional information. If none is given (default) marginal sums are used.}
}

\details{
The d' index is derived from Kulback-Leibler distance (as is Shannon's diversity index), and calculates how strongly a species deviates from a random sampling of interacting partners available. It ranges from 0 (no specialisation) to 1 (perfect specialist). In the case of a pollination web, a pollinator may be occurring only on one plant species, but if this species is the most dominant one, there is limited evidence for specialisation. Hence this pollinator would receive a low value. In contrast, a pollinator that occurs only on the two rarest plants would have a very high value of d'.

The idea of this index is laid out in Blüthgen et al. (2006). It basically calculates the Shannon-diversity for each column (delivering the raw d-values) and re-ranges them between the theoretical maximum and minimum (yielding values between 0 and 1). dmax is given analytically (see paper or code), but dmin must be found `heuristically', since the web can only contain integers. The idea behind the heuristic minimum is that d will be minimal when observed values differ least from expected values based on marginal distributions.

The way this function is implemented, it calculates expected values for each cell based on the product of observed marginal sums (i.e. column and row sums) times \code{sum(web)}. Then it rounds off to integers and allocates the remaining interactions in two steps: First, all columns and rows with marginal sums of 0 obtain one interaction into the cell with the highest expected value. Secondly, all remaining interactions are distributed according to difference between present and expected value: those cells with highest discrepancy receive an interaction until the sum of all entries in the new web equals those in the original web. Now the d-values for this web are calculated and used as dmin.

Simple rounding of expected values would lead to empty columns or rows, i.e. the dmin-web would be of lower dimension than the original web.

dfun returns the d' values for the lower trophic level. Use \code{fun(t(web))} to get the d'-values for the higher trophic level (as does 'specieslevel').

d' is one of several species-level network indices. It's generalisation to the entire interaction web is called H2'.

The \code{abundances} vector allows to incorporate independent estimates of the abundances of the HIGHER trophic level. In a pollination web, pollinator abundances may be very different from those estimated by the interaction matrix column sums. This has also, obviously, large consequences for the specialisation: A plant being pollinated by a bee that is common on this plant, but very rare in general, will show a low specialisation unless bee abundances are corrected for. Data given in the abundance vector are here used in replacement for the row sums, both in the d-function itself, as well as in the calculation of the minimum ds.
}

\value{
  \item{dprime}{d'-value for each species}
  \item{d}{Raw d-value for each species, i.e. before ranging it between 0 and 1.}
  \item{dmin}{Minimum d-value for each species, based on a perfect nesting of the matrix; see details.}
  \item{dmax}{Maximum d-value theoretically possible given the observed number of interactions and the observed marginal distributions.}
}

\references{ 
Blüthgen, N., Menzel, F. \& Blüthgen, N. 2006 Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 12
 }

\author{Carsten F. Dormann and Jochen Fründ}
\note{ 
As detailed above, deriving the dmin-values \sQuote{heuristically} leaves room to some variation. We are very happy with this implementation, but you may want to program something yourself ...
Furthermore, incorporating independent measures of abundance are straightforward for most of the code. However, I (CFD) do not know how to calculated the maximum d values in this case. So for the time being we have implemented a quick-fix.
}

\seealso{ \code{\link{H2fun}} for a similar function for the entire network. \code{\link{specieslevel}} for a method that, amongst other indices, calls \code{dfun}. }

\examples{
data(Safariland)
dfun(Safariland) # gives d-values for the lower trophic level
# now using independent abundance estimates for higher trophic level:
dfun(Safariland, abundances=runif(ncol(Safariland)))

dfun(t(Safariland)) #gives d-values for the higher trophic level
}

\keyword{ htest }

