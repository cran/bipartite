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
We only had three types of bipartite webs in mind when writing this package: seed-disperser, plant-pollinator and predator-prey systems. In how far it makes sense to use these
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
Version: \tab 0.8\cr
Date: \tab 2008-12-19\cr
License: \tab GPL \cr
}
}

\section{versionlog}{

0.8 (release date: 18-Dec-2008)
\item{bug fix in \code{\link{shuffle.web}}}{Function didn't work correctly when the web contained more rows than columns. (Thanks to Anna Traveset for spotting and reporting!)}
\item{minor modifications in networklevel}{Sometimes additional information is available and species can be included in a network, although they have no interactions with other species. In this case, one might want to use these species in the network, too. A new option (empty.web), allows to keep empty rows and columns, although for some functions an emptied web must be used (e.g. degree distributions). }
\item{bug fix in \code{\link{H2fun}}}{As in the last \code{\link{H2fun}} bug fix, sometimes H2 became negative.}
\item{bug (?) fix in \code{\link{networklevel}}}{Shannon diversity is based on the log of interactions. If this value is 0, as it is for most network entries, an NA is produced. As a consequence, Shannon's H (now also given in output) is based not on all interactions, but only those > 0. In consequence, Shannon evenness should also only be SH divided by the number of realised interactions (\code{log(sum(web>0))}). That is now the case.}

0.74 (release date: 24-Oct-2008)
\item{functional specialisation (\code{\link{functspec}}) bug fix}{Paths were double the true length, hence minimum was 2, rather than 1.}
\item{\code{\link{H2fun}} bug fix}{Since the search for H2min is heuristic, H2uncorr can sometimes be lower than H2min; in that case, H2fun returned a value greater 1, while it should be one exactly.}
\item{new function \code{\link{nestedness.corso}} }{Calculates (weighted) nestedness according to Corso et al. (2008) and Galeano et al. (2008).}
\item{new function \code{\link{discrepancy}}}{Calculates discrepancy according to Brualdi & Sanderson (1999), deemed to be best ever measure of nestedness; also gives an example for a binary null model analysis based on \pkg{vegan}'s \code{commsimulator}.}
\item{other}{Correction of several minor typos on the help pages; removal of "~" in help files; same citation style throughout; new cross references (especially for the nestedness functions); in \code{\link{networklevel}}, nestedness is now calculated using \pkg{vegan}'s \code{nestedtemp} due to matrix inversion problems reported for binmatnest.}

0.73 (release date: 1-Sept-2008)
\item{new feature}{\code{\link{plotweb}}}{Named abundance-vector for each level can be used.}
\item{new function \code{\link{plotweb2}} (not debugged!)}{For plotting tripartite networks.}

0.72 (release date: 12-June-2008)
\item{new function: functional specialisation with \code{\link{functspec}}}{See Dalgaard et al. (2008).}
\item{new function: interface to \pkg{sna} through \code{\link{as.one.mode}}}{Allows calculation of path lengths, centrality, betweenness and other indices developed for one-mode networks.}
\item{bug-fix}{Error in \code{\link{plotweb}} when no species labels were given.}

}

\author{
Carsten F. Dormann and Bernd Gruber, with additional code from Jochen Fründ

Maintainer: Carsten Dormann \email{carsten.dormann@ufz.de} and Bernd Gruber \email{bernd.gruber@ufz.de}

}
\references{
Bascompte, J., Jordano, P. and Olesen, J. M. (2006) Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. (2002) Quantitative descriptors of food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N., Menzel, F. and Blüthgen, N. (2006) Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 12

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen, N. (2007) Specialization, constraints, and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Corso G., de Araújo A.I.L. and de Almeida A.M. (2008) A new nestedness estimator in community networks. \emph{arXiv}, 0803.0007v1 [physics.bio-ph]

Dalsgaard, B., A. M. Martín González, J. M. Olesen, A. Timmermann, L. H. Andersen, and J. Ollerton. (2008) Pollination networks and functional specialization: a test using Lesser Antillean plant-hummingbird assemblages. \emph{Oikos} \bold{117}, 789--793

Galeano J., Pastor J.M. and Iriondo J.M. (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Memmott, J., Waser, N. M. and Price, M. V. (2004) Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Vazquez, D. P. and Aizen, M. A. (2004) Asymmetric specialization: A pervasive feature of plant-pollinator interactions. \emph{Ecology} \bold{85}, 1251--1257


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
