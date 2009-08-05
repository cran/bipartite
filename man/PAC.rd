\encoding{latin1}
\name{PAC}
\alias{PAC}

\title{Potential for Apparent Competition}
\description{
Quantifies, for each pair of lower trophic level species, the potential for showing apparent competition with another species, mediated through the higher trophic level.
}
\usage{
PAC(web)
}

\arguments{
  \item{web}{A host-parasitoid network (or alike), where the entries represent the \bold{sum of parasitoids emerging from the interactions between parasitoid and host} (i.e. number of interactions * number of parasitoid individuals emerging from each host). Only if there is only one parasitoid per host this web will be the same as that used in all other calculations in this package!}
}

\details{
Calculates the potential for apparent competition (Holt 1977), following the formula given in M�ller et al. (1999) and Morris et al. (2005). See also Morris et al. (2004) for an experimental test.
}

\value{
 Returns a k x k matrix with entries d.ij, where k is the number of species in the lower trophic level and i and j are lower trophic level species. Entries in the upper triangle represent the effect of the species j on species i (d.ij), i.e. the effect of the species in the column onto that in the row. In analogy, the lower triangle represent the inverse direction (d.ji), i.e. the effect of the row species onto the column species. Diagonal entries are \dQuote{apparent intraspecific competition}.
}

\note{
The idea is that in host-parasitoid networks one host also affects other hosts by the number of parasitoid that hatch from it and are thus added to the pool of parasitoids. An abundant, large host can (involuntarily) contribute many parasitoids to the pool, thus also increasing the parasitoid burden of other hosts. This looks like competition between the two hosts, while in fact it is mediated through the other trophic level. 

Whether this concept can be usefully applied to mutualist networks (such as flower visitation networks, aka pollination webs) is still under debate. The example below has thus to be seen as a technical, not a biological example.
}

\references{ 
Holt, R. D. 1977 Predation, apparent competition and the structure of prey communities. \emph{Theoretical Population Biology} \bold{12}, 197--229.

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2004 Experimental evidence for apparent competition in a tropical forest food web. \emph{Nature} \bold{428}, 310--313.

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2005 Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449--462.

M�ller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de}}

\seealso{
}

\examples{
data(Safariland)
PAC(Safariland)
}

\keyword{ package }

