\name{as.one.mode}

\alias{as.one.mode}

\title{Conversion of a network matrix}

\description{
 This helper function converts a bipartite matrix into a one-mode matrix.
}

\usage{
as.one.mode(web, fill = 0)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species 
    as columns and number of interactions as entries.}
  \item{fill}{What shall unobserved combinations be represented as in the one-mode matrix (see below)? Defaults to 0.}
}

\details{
  In bipartite (or: two-mode) networks, participants are of different types (e.g. pollinators and plants, 
  actors and parties in social research). Hence, a party cannot connect to another party except through 
  actors. A pollinator interacts with another pollinator only through the host plant.
  
  Much network theory, however, is based on one-mode networks, where all participants are listed in one vector,
  i.e. plants and pollinators alike, actors together with events. This function here transforms the more condensed
  bipartite representation into a one-mode-representation, filling the unobserved type of interactions (
  i.e. plants with plants and pollinators with pollinators) with 0 (unless you specify it differently in \option{fill}).

  The lower trophic level (e.g. plants or rows) is listed first, then the higher trophic level (e.g. pollinators or columns). 
  Hence, pollinator 2 becomes species number r+2, where r is the number of rows of the network matrix.
  
  The benefit of this conversion is access to the wonderful R-package Social Network Analysis (\pkg{sna}), 
  with its many one-mode indices (such as \code{\link[sna]{betweenness}}, \code{\link[sna]{closeness}}, \code{\link[sna]{centralization}}, 
  \code{degree}, \code{kpath.census} and so forth). Furthermore, \code{\link[sna]{gplot}} in that package also provides cool network depictions well worth checking out.
  
  With respect to \pkg{bipartite}, \code{as.one.mode} is employed in the function \code{\link{nodespec}}, which itself uses the \pkg{sna}-function \code{\link[sna]{geodist}}.
}

\value{
  A matrix of dimension (n+k) x (n+k), where n and k are the dimensions of the input 
  web. Both dimensions are given the names of the orignal web (first the lower, then
  the higher trophic level).
}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de} }

\note{ 
This function is not at all clever or complex. In fact, most users may not consider it worth 
providing at all.
}

\examples{
data(Safariland)
image(Safariland)
image(as.one.mode(Safariland))
}

\keyword{package}

