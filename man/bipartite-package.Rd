\encoding{latin1}
\name{bipartite-package}
\alias{bipartite-package}
\alias{bipartite}
\docType{package}
\title{
Analysis of bipartite ecological webs.
}
\description{
Bipartite provides functions to viualise webs and calculate a series of indices commonly
used to describe pattern in ecological webs. It focusses on webs consisting of only
two trophic levels, e.g. pollination webs or predator-prey-webs. Visualisation is
important to get an idea of what we are actually looking at, while the indices summarise
different aspects of the webs topology.}
\details{

We only had three types of bipartite webs in mind when writing this package: seed-disperser,
plant-pollinator and predator-prey systems. In how far it makes sense to use these
functionalities for other systems (or indeed for these systems) lies in the hands
of the user. Please refer to the literature cited for details on the theory behind
the indices.

Input for most analyses is an interaction matrix of m higher level species with n
lower level species, i.e. an n x m matrix, where higher trophic level species are
in columns, lower level in rows. Column and row names can be provided. This is fundamentally
different from multi-dimensional webs, which are organised as k x k matrix, i.e. each
species against each other. Such a format is incompatible with the functions we provide
here.

The first step is to visualise the interaction web. Two functions are on offer here:
one (\code{\link{visweb}}) simply plots the matrix in colours depicting the strength
of an interaction and options for re-arranging columns and rows (e.g. to identify
compartments or nesting). The other function (\code{\link{plotweb}}) plots the actual
web with participants (as two rows of rectangles) connected by lines (proportional to
interaction strength).

The second step is to calculate various indices describing network topography. There
are two different levels this can be achieved at: the entire web (using function
\code{\link{networklevel}}) or the individual species (using function \code{\link{specieslevel}}).
All other functions in the package are helpers, although some can be called on their
own and return the respective result (\code{\link{dfun}}, \code{\link{H2fun}} and
\code{\link{second.extinct}} with \code{\link{slope.bipartite}}).

See function description for details and examples.

\tabular{ll}{
Package: \tab bipartite\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2008-01-14\cr
License: \tab GPL \cr
}

}
\author{
Carsten F. Dormann and Bernd Gruber, with additional code from Jochen Fründ and Denis Lippok

Maintainer: Carsten Dormann \email{carsten.dormann@ufz.de} or Bernd Gruber \email{bernd.gruber@ufz.de}

}
\references{
Bascompte, J., Jordano, P. and Olesen, J. M. 2006 Asymmetric coevolutionary networks
facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. 2002 Quantitative descriptors of
food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N., Menzel, F. and Blüthgen, N. 2006 Measuring specialization in species
interaction networks. \emph{BMC Ecology} \bold{6}, 12

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen, N. 2007 Specialization,
constraints, and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Memmott, J., Waser, N. M. and Price, M. V. (2004) Tolerance of pollination networks to
species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. 2007 Habitat modification alters the
structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Vazquez, D. P. and Aizen, M. A. 2004 Asymmetric specialization: A pervasive feature of
plant-pollinator interactions. \emph{Ecology} \bold{85}, 1251--1257


}
\keyword{ package }
\seealso{
}
\examples{
data(Safariland)
plotweb(Safariland)
visweb(Safariland)
networklevel(Safariland)
specieslevel(Safariland)
}
