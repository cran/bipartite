\encoding{latin1}
\name{networklevel}
\alias{networklevel}

\title{Analysis of bipartite webs at the level of the entire network}
\description{
Calculates a variety of indices and values for a bipartite network.
}
\usage{
networklevel(web, index="ALL", ISAmethod="Bluethgen", SAmethod = "Bluethgen",
    extinctmethod = "r", nrep = 100, plot.it.extinction = FALSE, plot.it.dd=FALSE, CCfun=median, dist="horn", normalise=TRUE, nest.weighted=FALSE, empty.web=TRUE)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{index}{One or more of the following (exact match only!):\dQuote{number of species},
  \dQuote{links per species}, \dQuote{connectance}, \dQuote{linkage density},\dQuote{web asymmetry},
  \dQuote{number of compartments}, \dQuote{generality},\dQuote{vulnerability},
  \dQuote{interaction evenness},\dQuote{compartment diversity}, \dQuote{cluster coefficient},
  \dQuote{H2},\dQuote{ISA},\dQuote{SA},\dQuote{extinction slope}, \dQuote{degreedistribution},
  \dQuote{niche overlap}, \dQuote{mean number of shared hosts}, \dQuote{togetherness},
  \dQuote{C-score}, \dQuote{V-ratio}, \dQuote{nestedness}, \dQuote{nestedness.corso},\dQuote{discrepancy}. \dQuote{ALL} (default) calculates all indices.}
  \item{ISAmethod}{Method to use for calculating interaction strength (= dependence)
  asymmetry; original by \dQuote{Bascompte}  is biased for singletons and few interactions
  (range 0 to infty); \dQuote{Bluethgen} (default) excludes singletons and corrects
  for low number of interactions (range -1 to 1); see Bascompte et al. (2006) and
  Blüthgen et al. (2007).}
  \item{SAmethod}{Which method to use to calculate the specification asymmetry: mean
  of log-transformed dependencies (\dQuote{log}) or Blüthgen's abundance-weighted mean
  (default); see Blüthgen et al. (2007).}
  \item{extinctmethod}{Specifies how species are removed from matrix: \dQuote{random}
  or \dQuote{abundance} (partial matching), where abundance removes species in the
  order of increasing abundance (i.e. rarest first); see Memmott et al. (2004).}
  \item{nrep}{Number of replicates for the extinction sequence analysis.}
  \item{plot.it.extinction}{logical; plot the extinction sequence graph; defaults to FALSE. }
  \item{plot.it.dd}{logical; plot the degree distribution fits?; defaults to FALSE. }
  \item{CCfun}{Method to use when calculating the clustering coefficient. Originally proposed as mean of cluster coefficients for each node. Defaults to \dQuote{median}, because cluster coefficients are strongly skewed. }
  \item{dist}{Distance metric to be used to calculate niche overlap; defaults to Horn's index, which is the recommendation of Krebs (Ecological Methodology); for other options see \code{vegdist} in \pkg{vegan}.}
  \item{normalise}{Logical; shall the C-score and togetherness metrics be normalised to a range of 0 to 1? Defaults to \code{TRUE}.}
  \item{nest.weighted}{Shall the interaction-weighted nestedness be calculated (see \code{\link{nestedness.corso}} for details)? Defaults to \code{FALSE}.}
  \item{empty.web}{Shall the empty columns and rows be deleted? Defaults to TRUE.}
}

\details{
This function implements a variety of the many (and still procreating) indices describing network topography. Some are embaracingly simple (such as number of species in each trophic level or the number of links (= non-zero cells) in the web). Others are variations on Shannon's diversity index applied to within column or within rows. Only extinction slope is newly implemented here, and hence described in a bit more detail.

Extinction slope works on a repeated random sequence of species extinctions (within one trophic level), and calculates the number of secondary extinctions (in the other level). These values are then averaged (over the \option{nrep} runs) and plotted against the number of species exterminated. The proportion still recent (on the y-axis) regressed against the proportion exterminated (on the x-axis) is hence standardised to values between 0 and 1 each. Through this plot, a hyperbolic regression is fitted, and the slope of this regression line is returned as an index of extinction sensitivity. The larger the slope, the later the extinction takes its toll on the other trophic level, and hence the higher the redundancy in the trophic level under consideration. Using \option{plot.it=F} also returns the graphs (set history to recording in the plotting window). Changing the \option{extinctionmethod} to \dQuote{abundance} will always result in the same sequence (by increasing abundance) and hence does not require replication.

Most indices are straightforward, one-line formulae; some, such as H2', also require a re-arranging of the matrix; and one, secondary extinction slope, internally requires iterative runs, making the function relatively slow. If you are not interested in the secondary extinction slopes, simply set \option{nrep=1} to make it much faster.
}

\value{
Depending on the selected indices, some or all of the below:
  \item{number of higher trophic species}{Number of species in the higher trophic level, i.e. \code{ncol(web)}.}
  \item{number of lower trophic species }{Number of species in the lower trophic level, i.e. \code{nrow(web)}.}
  \item{links per species}{Mean number of links per species (qualitative): sum of links divided by number of species.}
  \item{connectance}{Realised proportion of possible links (Dunne et al. 2002): sum of links divided by number of cells in the matrix (= number of higher times number of lower trophic level species). This is the standardised 
  number of species combinations often used in co-occurrence analyses (Gotelli \& Graves 1996)}
  \item{linkage density}{Mean number of interactions per species (quantitative); see Tylianakis et al. (2007).}
  \item{web asymmetry}{Balance between numbers in the two levels: positive numbers indicate more lower-trophic level species, negative more higher-trophic level species; implemented as \code{(nrow(web)-ncol(web))/sum(dim(web))}; web asymmetry is a null-model for what one might expect in dependence asymmetry: see Blüthgen et al. (2007).}
  \item{number of compartments}{Compartments are sub-sets of the web which are not connected (through either higher or lower trophic level) to another compartment. Mathematically, they are Jordan blocks, but this implementation is rule-based (and fast). They are also nicely visualised in the \code{\link{visweb}} function.}
  \item{generality}{Effective mean number of prey per predator; see Tylianakis et al. (2007).}
  \item{vulnerability}{Effective mean number of predator per prey; see Tylianakis et al. (2007).}
  \item{Shannon diversity}{Shannon's diversity of interactions (i.e. network entries).}
  \item{interaction evenness}{Shannon's evenness for the web entries, treating zeros as no data.}
  \item{Alatalo interaction evenness}{A different measure for web entry evenness, as proposed by Müller et al. (1999). }
  \item{compartment diversity C.D.}{Shannon's diversity index across compartment sizes (i.e. number of participants); see Tylianakis et al. (2007).}
  \item{cluster coefficient}{The CC for a network is the average CC of its members. CC for each node, in turn, is simply the number of realised links devided by the number of possible links. Introduced by Watts \& Strogatz (1998) and described in Wikipedia under \url{http://en.wikipedia.org/w/index.php?title=Clustering_coefficient}.}
  \item{H2}{H2' is a network-level measure of specialisation. It ranges between 0 (no specialisation) and 1 (complete specialisation). H2' is a measure of discrimination, i.e. calculated in comparison of no specialisation (see \code{\link{H2fun}} for details. To avoid confusion of keys (apostrophe vs. accent), we call the H2' only H2 here.}
  \item{dependence asymmetry}{Explaining dependence asymmetry will require more space than just a few lines. In essence, it is also a measure of specialisation, across both trophic levels. Proposed by Bascompte et al. (2006) and critiqued and alterations proposed by Blüthgen et al. (2007). The latter also show that dependence asymmetry can be almost entirely explained by web asymmetry (see above). Positive values (only possible of \option{DAmethod="Bluethgen"}) indicate higher dependence in the higher trophic level.}
  \item{specialisation asymmetry}{Asymmetry (higher vs. lower trophic level) of specialisation now based on d' (see \code{\link{dfun}}), which is insensitive to the dimensions of the web. Again, two options of calculation are available: the one proposed by Blüthgen et al. (2007), where they weight the specialisation value for each species by its abundance (\option{SAmethod="Bluethgen"}) or where d'-values are log-transformed (argueing that d'-values are indeed log-normally distributed: \option{SAmethod="log"}). Since the mean d-value for the lower trophic level is subtracted from that of the higher, positive values indicate a higher specialisation of the higher trophic level.}
  \item{extinction slope higher trophic level}{Slope of the secondary extinction sequence in the higher trophic level, following extermination of species in the lower trophic level; see details.}
  \item{extinction slope lower trophic level}{Slope of the secondary extinction sequence in the lower trophic level, following extermination of species in the higher trophic level; see details.}
  \item{degree distribution}{Coefficients and fits for three different functions to degree distributions: exponential, power law and truncated power law. See \code{\link{degreedistr}} for details and references.}
  \item{niche overlap}{Mean similarity in interaction pattern between species of the same trophic level, calculated by default as Horn's index. Values near 0 indicate no common use of niches, 1 indicates perfect niche overlap. (In how far it makes sense for hosts of predators to call their commonness in enemies \dQuote{niche overlap} is a different issue. There are people calling predators \dQuote{negative resources} (couldn't be asked to look up the reference). I would probably rather call it similarity in predators.)}
  \item{mean number of shared hosts}{The simplest measure of co-occurrence and hence similarity in host preferences; based on Roberts \& Stone (1990) and Stone \& Roberts (1992).}
  \item{togetherness}{Mean number of co-occupancies across all species-host-combinations; the whole matrix is scanned for submatrices of the form c(0,0,1,1), representing perfect matches of co-presences and co-absences. These are counted for each pairwise species combination, and averaged. Based on Stone \& Roberts (1992).}
  \item{C-score}{Mean (normalised) number of checkerboard combinations across all higher trophic level species. Values close to 1 indicate that there is evidence for disaggregation, e.g.~through competition. Value close to 0 indicate aggregation of species (i.e.~no repelling forces between species). See Stone and Roberts (1990) for details.}
  \item{V-ratio}{Variance-ratio of species numbers to individual numbers within species for the higher trophic level. Values larger than 1 indicate positive aggregation, values between 0 and 1 indicate disaggregation of species. See Schluter (1984) for details.}
  \item{nestedness}{Nestedness-temperature of the matrix. For details see \code{\link{nestedness}} and Rodríguez-Gironés \& Santamaria (2002). Notice that the function \code{\link{nestedness}}, as called by \code{networklevel}, does not calculate any null model, simply because it is too computer-intensive. If you are interested in the different null models, please use the function \code{\link{nestedness}} or \code{nestedtemp} in \pkg{vegan} directly.}
  \item{nestedness.corso}{Another way to calculate nestedness, proposed by Corso et al. (2007). A version that considers interaction frequencies (and is hence weighted) has been proposed (Galeano et al. 2007), and its value is returned for option \option{nest.weighted=TRUE}. Both range between 0 (perfect nestedness) and 1 (perfect chaos) and have not been properly evaluated!}
  \item{discrepancy}{Discrepancy as proposed by Brualdi \& Sanderson (1999); see \code{[\link{discrepancy}} for details.}
}

\note{Some indices require rather long computation times on large webs. Particularly slow (for large webs) are: \dQuote{extinction slope} and \dQuote{nestedness}. For ordinary-sized webs (tens by tens), \code{networklevel} is reasonably fast. An example of a large web is \code{kato1990}: try it, possibly using \code{system.time} as a wrapper, and prepare to do some reading on the side ...}

\references{ 
Bascompte, J., Jordano, P. and Olesen, J. M. 2006 Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen N. 2007 Specialization, constraints and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Corso G, de Araújo AIL, de Almeida AM (2008) A new nestedness estimator in community networks. \emph{arXiv} 0803.0007v1 [physics.bio-ph]

Dunne, J. A., R. J. Williams, and N. D. Martinez. 2002 Food-web structure and network theory: the role of connectance and size. \emph{Proceedings of the National Academy of Science USA} \bold{99}, 12917--12922

Galeano J, Pastor JM, Iriondo JM (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Gotelli, N. J., and G. R. Graves. 1996 Null Models in Ecology. Smithsonian Institution Press, Washington D.C.

Memmott, J., Waser, N. M. and Price M. V. 2004 Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Roberts, A. and Stone, L. 1990 Island-sharing by archipelago species. \emph{Oecologia} \bold{83}, 560--567

Rodríguez-Gironés M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935

Schluter, D. 1984 A variance test for detecting species associations, with some example
applications. \emph{Ecology} \bold{65}, 998-1005.

Stone, L. and Roberts, A. 1990 The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79.

Stone, L. and Roberts, A. 1992 Competitive exclusion, or species aggregation? An aid in deciding. \emph{Oecologia} \bold{91}, 419--424

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. 2007 Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Watts, D. J. and Strogatz, S. 1998 Collective dynamics of \sQuote{small-world} networks. \emph{Nature} \bold{393}, 440--442
}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de}}

\seealso{ Some functions are implemented separately: \code{\link{H2fun}}, \code{\link{second.extinct}}, \code{\link{degreedistr}}, \code{\link{C.score}} and \code{\link{V.ratio}}
}

\examples{
data(Safariland)
networklevel(Safariland)
}

\keyword{ package }


