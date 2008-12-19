\encoding{latin1}

\name{nestedness.corso}

\alias{nestedness.corso}

\title{Calculates the (weighted) nestedness of a network according to the method proposed by Corso et al. (2008), with weighting option by Galeano et al. (2008).}

\description{Calculates a nestedness value (not unlike the matrix temperature returned by \code{\link{nestedness}}, but in the interval [0,1] (0 means no nesting, 1 perfect nesting). Nestedness according to this function differs fundamentally in the way it is calculated. Roughly, it describes the observed distances between links in the network matrix, and compares that to what is possible under the given connectance.

 The weighted version builds on the fact that at some point during the calculations of Corso et al.'s nestedness index the Manhattan distance is used to calculate distances between cells. Here, we can slot in a simple weight (in this case the respective dependencies) and hence weight the distance by the number of observations (or rather the dependencies).
}

\usage{
nestedness.corso(web, weighted=FALSE, reps=5000)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species
    as columns and number of interactions as entries.}
  \item{weighted}{Should the quantitative version, using number of observations to calculate weights, be used? Defaults to \code{FALSE}.}
  \item{reps}{Number of replicates for constructing random networks (used in the weighted version to find the minimum nestedness). Defaults to 5000. Anything less than 1000 would be very unsensible, while 10000 and more will take up time (only a few seconds per 1000).}
}

\details{
 The logic and descriptions are available from the two references. However, it may be useful to point out that while there is a straightforward way to calculate observed, maximum and minimum nestedness for the unweighted version, the weighted version has an element of heuristics (i.e. fiddling).
 
 In principle, both versions sort the matrix by arranging rows and columns to densest packing, then calculate distances between (occupied) cells and sum those. Then the same is done for the minimally and maximally nested matrix with the same dimensions and number of links. For the weighted version, however, only the maximally nested can be identified uniquely, not the unnested maximum. Hence, we resorted to using the \code{\link{shuffle.web}} function to produce random, high-entropy webs, that serve to delimit the lowest possible nestedness.
}

\value{
  The nestedness value, ranging between 0 (chaos) and 1 (perfect nesting).
}

\references{ 

Corso G, de Araújo AIL, de Almeida AM (2008) A new nestedness estimator in community networks. \emph{arXiv} 0803.0007v1 [physics.bio-ph]

Galeano J, Pastor JM, Iriondo JM (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Ulrich, W. and Gotelli, N.J. (2007) Null model analysis of species nestedness patterns. \emph{Ecology} \bold{88}, 1824--1831


 }

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de}}

\note{ While the logic behind the unweighted index seems at least plausible, I am not so sure about the weighted version. As far as I know, there is no systematic comparison with other nestedness indices, nor has it been shown that the indices are \dQuote{well-behaved}. Thus, either the esteemed user writes a paper along these lines, or he/she uses the indices with caution (i.e. keep the findings to themselves, because nobody would know what to make of them).

 For a well-performed comparison of other nestedness indices see Ulrich \& Gotelli (2007).

}

\seealso{ See also as \code{\link{nestedness}} and \code{\link{discrepancy}}}

\examples{
data(Safariland)
nestedness.corso(Safariland, weighted=TRUE, reps=1000)
}

\keyword{package}
