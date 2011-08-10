\encoding{UTF-8}
\name{PDI}
\alias{PDI}

\title{Paired Differences Index}
\description{ Computes the Paired Differences Index }

\usage{
PDI(web, log=FALSE)
}

\arguments{
  \item{web}{A bipartite interaction web, i.e.~a matrix with higher (cols) and lower (rows) trophic levels.}
  \item{log}{logical; since number of interactions is often highly skewed, the log yields a more even spread of PDI-values across species. Defaults to \code{FALSE}.}
}

\details{
There are many ways to skin a cat. This is a more recent addition proposed by Poisot et al. (2011a) and used for mutualistic network analysis by Poisot et al. (2011b). This function can be used alone or through \code{\link{specieslevel}}. 

If P is the proportion of interactions a species engages in, then PDI for this species is computed as:

 sum(P1 - Pi) /(H-1), 
 
where P1 is the highest number of interactions in a link, while Pi are the remaining values. H is the number of interactors (e.g. plant species if the target species is a pollinators). 

As so many other indices, also PDI is somewhat dubious for species which were observed only once (singletons). In this case P1 is 1 (all interactions are on one plant species), all Pi are 0 and hence PDI is 0, too. This is not logical. For a singleton we cannot tell its specialisation, only its discrimination (i.e.whether it happens to visit a common or rare plant species), as is done by d' (implemented in \code{\link{dfun}}).

}

\note{
This index was not included in the comparison of Dormann (2011).
}

\value{
Returns a vector with PDI values between 0 (generalist) and 1 (perfect specialist).
}

\references{
Dormann, C.F. (2011) How to be a specialist? Quantifying specialisation in pollination networks. \emph{Network Biology} \bold{1}, 1--20

Poisot, T., Lepennetier, G., Martinez, E., Ramsayer, J., and Hochberg, M.E. (2011a) Resource availability affects the structure of a natural bacteria-bacteriophage community. \emph{Biology Letters} \bold{7}, 201--204

Poisot, T., Bever, J.D., Nemri, A., Thrall, P.H., and Hochberg, M.E. (2011b) A conceptual framework for the evolution of ecological specialisation. \emph{Ecology Letters} \bold{14}, 841--851
}

\author{ Carsten F. Dormann <carsten.dormann@ufz.de> }

\seealso{ See also \code{\link{specieslevel}}. }

\examples{
data(Safariland)
PDI(Safariland) # for pollinators
PDI(t(Safariland), log=TRUE) # for plants
}
\keyword{ package}

