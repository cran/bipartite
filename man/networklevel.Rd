\encoding{latin1}
\name{networklevel}
\alias{networklevel}

\title{Analysis of bipartite webs at the level of the entire network}
\description{
Calculates a variety of indices and values for a bipartite network.
}
\usage{
networklevel(web, index="ALL", ISAmethod="Bluethgen", SAmethod = "Bluethgen",
    extinctmethod = "r", nrep = 100, plot.it.extinction = FALSE, plot.it.dd=FALSE,
    CCfun=median, dist="horn")
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{index}{One or more of the following (exact match only!):\dQuote{number of species},
  \dQuote{number of links},\dQuote{linkage density},\dQuote{web asymmetry},\dQuote{number of compartments},
  \dQuote{generality},\dQuote{vulnerability},\dQuote{interaction evenness},\dQuote{compartment diversity},
  \dQuote{cluster coefficient},\dQuote{H2},\dQuote{ISA},\dQuote{SA},\dQuote{extinction slope},
  \dQuote{niche overlap}. \dQuote{ALL} (default) calculates all indices.}
  \item{ISAmethod}{Method to use for calculating interaction strength (= dependence)
  asymmetry; original by \dQuote{Bascompte}  is biased for singletons and few interactions
  (range 0 to infty); \dQuote{Bluethgen} (default) excludes singletons and corrects
  for low number of interactions (range -1 to 1); see Bascompte et al. (2006) and
  Blüthgen et al. (2007).}
  \item{SAmethod}{Which method to use to calculate the specification asymmetry: mean
  of log-transformed dependences (\dQuote{log}) or Blüthgen's abunance-weighted mean
  (\dQuote{Bluethgen}: default); see Blüthgen et al. (2007).}
  \item{extinctmethod} {Specifies how species are removed from matrix: \dQuote{random}
  or \dQuote{abundance} (partial matching), where abundance removes species in the
  order of increasing abundance (i.e. rarest first); see Memmott et al. (2004).}
  \item{nrep}{Number of replicates for the extinction sequence analysis.}
  \item{plot.it.extinction}{logical; plot the extinction sequence graph; defaults to FALSE. }
  \item{plot.it.dd}{logical; plot the degree distribution fits?; defaults to FALSE. }
  \item{CCfun}{Method to use when calculating the clustering coefficient. Originally proposed as mean of cluster coefficients for each node. Defaults to \dQuote{median}, because cluster coefficients are strongly skewed. }
  \item{dist}{ Distance metric to be used to calculate niche overlap; defaults to Horn's index, which is the recommendation of Krebs (Ecological Methodology); for other options see \code{vegdist} in package \code{vegan}.}
}

\details{
This function implements a variety of the many (and still procreating) indices describing network topography. Some are embaracingly simple (such as number of species in each trophic level or the number of links (= non-zero cells) in the web). Others are variations on Shannon's diversity index applied to within column or within rows. Only extinction slope is newly implemented here, and hence described in a bit more detail.

Extinction slope works on a repeated random sequence of species extinctions (within one trophic level), and calculates the number of secondary extinctions (in the other level). These values are then averaged (over the \option{nrep} runs) and plotted against the number of species exterminated. The proportion still recent (on the y-axis) regressed against the proportion exterminated (on the x-axis) is hence standardised to values between 0 and 1 each. Through this plot, a hyperbolic regression is fitted, and the slope of this regression line is returned as an index of extinction sensitivity. The larger the slope, the later the extinction takes its toll on the other trophic level, and hence the higher the redundancy in the trophic level under consideration. Using \option{plot.it=F} also returns the graphs (set history to recording in the plotting window). Changing the \option{extinctionmethod} to \dQuote{abundance} will always result in the same sequence (by increasing abundance) and hence does not require replication.

Most indices are straightforward, one-line formulae; some, such as H2', also require a re-arranging of the matrix; and one, secondary extinction slope, internally requires iterative runs, making the function relatively slow. If you are not interested in the secondary extinction slopes, simply set \option{nrep=1} to make it much faster.
}

\value{
  \item{number of higher trophic species}{Number of species in the higher trophic level, i.e. \code{ncol(web)}.}
  \item{number of lower trophic species }{Number of species in the lower trophic level, i.e. \code{nrow(web)}.}
  \item{mean number of links}{Mean number of links per species (qualitative).}
  \item{linkage density}{Mean number of interactions per species (quantitative); see Tylianakis et al. (2007).}
  \item{web asymmetry}{Balance between numbers in the two levels: positive numbers indicate more lower-trophic level species, negative more higher-trophic level species; implemented as \code{(nrow(web)-ncol(web))/sum(dim(web))}; web asymmetry is a null-model for what one might expect in dependence asymmetry: see Blüthgen et al. (2007).}
  \item{number of compartments}{Compartments are sub-sets of the web which are not connected (through either higher or lower trophic level) to another compartment. Mathematically, they are Jordan blocks, but this implementation is rule-based (and fast). They are also nicely visualised in the \code{\link{visweb}} function.}
  \item{generality}{Mean number of prey per predator; see Tylianakis et al. (2007).}
  \item{vulnerability}{Mean number of predator per prey; see Tylianakis et al. (2007).}
  \item{interaction evenness}{Shannon's evenness for the web entries, treating zeros as no data.}
  \item{Alatalo interaction evenness}{A different measure for web entry evenness, as proposed by Müller et al. (1999). }
  \item{compartment diversity C.D.}{Shannon's diversity index across compartment sizes (i.e. number of participants); see Tylianakis et al. (2007).}
  \item{cluster coefficient}{The CC for a network is the average CC of its members. CC for each node, in turn, is simply the number of realised links devided by the number of possible links. Introduced by Watts \& Strogatz (1998) and described in Wikipedia under \url{http://en.wikipedia.org/w/index.php?title=Clustering_coefficient}.}
  \item{H2}{H2' as a network-level measure of specialisation. It ranges between 0 (no specialisation) and 1 (complete specialisation). H2' is a measure of discrimination, i.e. calculated in comparison of no specialisation (see \code{\link{H2fun}} for details. To avoid confusion of keys (apostrophe vs. accent), we call the H2' only H2 here.}
  \item{dependence asymmetry}{Explaining dependence asymmetry will require more space than just a few lines. In essence, it is also a measure of specialisation, across both trophic levels. Proposed by Bascompte et al. (2006) and critiqued and alterations proposed by Blüthgen et al. (2007). The latter also show that dependence asymmetry can be almost entirely explained by web asymmetry (see above). Positive values (only possible of DAmethod=\dQuote{Bluethgen}) indicate higher dependence in the higher trophic level.}
  \item{specialisation asymmetry}{Asymmetry (higher vs. lower trophic level) of specialisation now based on d' (see \code{\link{dfun}}), which is insensitive to the dimensions of the web. Again, two options of calculation are available: the one proposed by Blüthgen et al. (2007), where they weight the specialisation value for each species by its abundance (\code{SAmethod=\dQuote{Bluethgen}}) or where d'-values are log-transformed (argueing that d'-values are indeed log-normally distributed: \code{SAmethod=\dQuote{log}}). Since the mean d-value for thelower trophic level is subtracted from that of the higher, positive values indicate a higher specialisation of the higher trophic level.}
  \item{extinction slope higher trophic level}{Slope of the secondary extinction sequence in the higher trophic level, following extermination of species in the lower trophic level; see details.}
  \item{extinction slope lower trophic level}{Slope of the secondary extinction sequence in the lower trophic level, following extermination of species in the higher trophic level; see details.}
  \item{niche overlap}{Mean similarity in interaction pattern between species of the same trophic level, calculated by default as Horn's index. Values near 0 indicate no common use of niches, 1 indicates perfect niche overlap. (In how far it makes sense for hosts of predators to call their commonness in enemies \dQuote{niche overlap} is a different issue. There are people calling predators \dQuote{negative resources} (couldn't be asked to look up the reference). I would probably rather call it similarity in predators.)}
}

\references{ 
Bascompte, J., Jordano, P. and Olesen, J. M. 2006 Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen N. 2007 Specialization, constraints and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Memmott, J., Waser, N. M. and Price M. V. 2004 Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. 2007 Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Watts, D. J. and Strogatz, S. 1998 Collective dynamics of \sQuote{small-world} networks. \emph{Nature} \bold{393}, 440--442
}

\author{ Carsten F. Dormann }

\seealso{ Two functions are implemented separately: \code{\link{H2fun}} and \code{\link{second.extinct}}.}

\examples{
data(Safariland)
networklevel(Safariland)
}

\keyword{ package }
