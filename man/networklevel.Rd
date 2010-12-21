\encoding{UTF-8}
\name{networklevel}
\alias{networklevel}

\title{Analysis of bipartite webs at the level of the entire network}
\description{
	Calculates a variety of indices and values for a bipartite network
}
\usage{
networklevel(web, index="ALL", ISAmethod="Bluethgen", 
    SAmethod = "Bluethgen", extinctmethod = "r", nrep = 100, 
    plot.it.extinction = FALSE, plot.it.dd=FALSE, CCfun=median, 
    dist="horn", normalise=TRUE, empty.web=TRUE, logbase="e", 
    intereven="sum", H2_integer=TRUE)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). 
Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{index}{One or more of the following (exact match only!):
    \option{number of species}, \option{connectance}, \option{web asymmetry}, \option{links per species}, 
    \option{number of compartments}, \option{compartment diversity}, \option{cluster coefficient}, \option{weighted cluster coefficient}, \option{degree 
distribution}, \option{mean number of shared partners}, \option{togetherness}, \option{C score}, \option{V ratio}, \option{discrepancy}, \option{nestedness}, 
\option{weighted nestedness}, \option{weighted NODF}, \option{extinction slope}, \option{robustness},
    \option{ISA},\option{SA}, \option{niche overlap},
    \option{generality},\option{vulnerability}, \option{linkage density}, \option{Fisher alpha}, \option{mean interaction diversity}, \option{interaction evenness}, 
\option{Alatalo interaction evenness}, \option{diversity}, \option{H2};
    \option{ALL} (default) calculates all indices (returning a list (since degree distribution fits are a table within this list and cannot be simplified)), 
    \option{ALLBUTDD} calculates all indices except degree distribution fits. This latter has the advantage that the output can be simplified to a vector; 
    \option{info} returns more general information on the network;
    \option{binary} returns a best-of selection of indices based on a binary network;
    \option{quantitative} returns a best-of selection of indices based on quantitative networks;
    \option{topology} returns indices more abstractly describing network properties. Also CHECK details below!
    }
  \item{ISAmethod}{Method to use for calculating interaction strength (= dependence) asymmetry; original by \option{Bascompte} is biased for singletons and 
few interactions (range 0 to infty); \option{Bluethgen} (default) excludes singletons and corrects for low number of interactions (range -1 to 1); see Bascompte 
et al. (2006) and Bl¸thgen et al. (2007).}
  \item{SAmethod}{Which method to use to calculate the specification asymmetry: mean of log-transformed dependencies (\option{log}) or Blüthgen's abundance-weighted 
mean (default); see Blüthgen et al. (2007).}
  \item{extinctmethod}{Specifies how species are removed from matrix: \option{random} or \option{abundance} (partial matching), where abundance removes species in the 
order of increasing abundance (i.e. rarest first); see Memmott et al. (2004).}
  \item{nrep}{Number of replicates for the extinction sequence analysis.}
  \item{plot.it.extinction}{logical; plot the extinction sequence graph; defaults to FALSE. }
  \item{plot.it.dd}{logical; plot the degree distribution fits?; defaults to FALSE. }
  \item{CCfun}{Method to use when calculating the clustering coefficient. Originally proposed as mean of cluster coefficients for each node. Defaults to \option{median}, 
because cluster coefficients are strongly skewed. }
  \item{dist}{Distance metric to be used to calculate niche overlap; defaults to Horn's index, which is the recommendation of Krebs (Ecological Methodology); for other options 
see \code{vegdist} in \pkg{vegan}.}
  \item{normalise}{Logical; shall the C-score and togetherness metrics be normalised to a range of 0 to 1? Defaults to \code{TRUE}.}
  \item{empty.web}{Shall the empty columns and rows be deleted? Defaults to TRUE.}
  \item{logbase}{Shall the various diversity indices be calculated to the base of e (default) or 2? Log2 is the proposal for generality and vulnerability by Bersier et al. (2002), 
while Shannon uses ln. The choice of the base will not affect the results qualitatively, only by a scaling factor.}
  \item{intereven}{Shall all cells of the matrix be used to calculate the interaction evenness (\option{intereven="prod"})? Or, as given by Bersier et al. (2002) and Tylianakis 
et al. (2007), should only the realised links be used (\option{intereven="sum"}; default)? Prod and sum refer to using the log of the \bold{product} of matrix dimensions (i.e. all cells) or the log of the \bold{sum} of non-zero cells (i.e. number of links)  as denominator in the evenness formula. See last paragraph of the details-section for views on these two options!}
  \item{H2_integer}{Logical; indicates whether values in web are integers. Passed on to \code{\link{H2fun}}; see there for details.}
}

\details{
This function implements a variety of the many (and still procreating) indices describing network topography. Some are embarrassingly simple and mere descriptors of a 
network's outer appearance (such as number of species in each trophic level or the number of links (= non-zero cells) in the web). Others are variations on Shannon's 
diversity index applied to within column or within rows. Only extinction slope is newly implemented here, and hence described in a bit more detail.

Note that Bersier et al. (2002) have three levels of values for some of their indices: qualitative (i.e. based on binary networks), quantitative (based on networks with 
information on the number of interactions observed for each link), and weighted-quantitative (where also independent information on the abundances of the participants 
is fed into the formula). At present, we implement only qualitative and quantitative indices! One exception is d', which can be called with an abundance vector (in function 
\code{\link{dfun}}). You can also get the qualitative version of quantitative indices (such as vulnerability) by simply calling networklevel on ``binarised'' data: 
\code{networklevel(Safariland>0, index="vulnerability")}. (Why you would want to do that, however, is currently beyond me.)(Thanks to Jason Tylianakis for proposing to 
put this clarification into the help!)

Integers or continuous values - what are the quantities in quantitative webs? Some web metrics expect in their typical formulation that the entries in the web-matrix are integers 
- e. g. H2' is defined relative to minimum and maximum based on marginal totals. Standard calculation (as in Bl¸thgen et al. 2006) uses an algorithm assuming values can only 
be integers. If your quantities are not constrained to be integers, the results may be wrong. Multiplication and rounding may or may not give consistent results, depending on
 rounding errors and the factor applied. Multiplication with high numbers such as 10 000 seem to be OK. For H2' a simplified calculation applicable to continuous numbers is 
available (by declaring option \option{H2_integer=FALSE} in \code{\link{H2fun}}). Note that values of H2' based on integers are not directly comparable to H2' based on continuous 
values (for sparse webs, H2'_continuous is much higher than H2'_integer). We are not sure how other indices are affected by non-integer values or by multiplication and rounding! 
Please let us know your experience.

Extinction slope works on a repeated random sequence of species extinctions (within one trophic level), and calculates the number of secondary extinctions (in the other level). 
These values are then averaged (over the \option{nrep} runs) and plotted against the number of species exterminated. The proportion still recent (on the y-axis) regressed against 
the proportion exterminated (on the x-axis) is hence standardised to values between 0 and 1 each. Through this plot, a hyperbolic regression is fitted, and the slope of this 
regression line is returned as an index of extinction sensitivity. The larger the slope, the later the extinction takes its toll on the other trophic level, and hence the higher the 
redundancy in the trophic level under consideration. Using \option{plot.it=F} also returns the graphs (set history to recording in the plotting window). Changing the 
\option{extinctionmethod} to \dQuote{abundance} will always result in the same sequence (by increasing abundance) and hence does not require replication.

Most indices are straightforward, one-line formulae; some, such as H2', also require a re-arranging of the matrix; and one, secondary extinction slope, internally requires
 iterative runs, making the function relatively slow. If you are not interested in the secondary extinction slopes, simply set \option{nrep=1} to make it much faster.

\bold{Sum or Prod: How to calculate interaction evenness?} I shall first put down my argument for ``prod'' and then Jason Tylianakis' arguments for ``sum''. 

Carsten: ``I do not want to defend a position I cannot hold against the flood of qualified criticism, and shall be happy to change the default to option ``sum'' (i.e. Jason's proposal). 
Nevertheless, I shall make a very brief attempt to defend my (and Nico's point of view). Imagine a completely different situation: I have ``counted'' birds in a landscape. From a more 
meticulous colleague I know that there are 27 bird species breeding at the moment, but on that two mornings that I went out, I could only hear 15. Now I want to calculate the Shannon 
diversity (and evenness) of birds in that landscape. The ``normal'' (in the sense of established) approach to use the data from my 15 species. But hold on: I KNOW there are more species 
out there. I don't know how many (i.e. there may be more than the 27 my colleague has found), but there are at least 27. If I only use the data from my 15 species, I will get a higher 
evenness value than when I also include the 12 zeros. My conclusion would be: I don't want to overestimate evenness only because I couldn't look long enough, thus I use all 27 values.''

Jason: ``I would disagree because what you ``know'' is based on your meticulous colleague's `sampling', which will also have its limits. If all you wanted was to know the total number 
of species there (assuming none have gone extinct), then what you propose is fine. However, the problem comes when you want to compare sites, and then sampling effort should be 
standardised. In most cases we know we don't have a full representation of the diversity (or food web) of an area, but we know for a given spatial or temporal sampling scale that one 
site differs from another in certain ways, and to me that is the most important. Anyway, it is all a question of scale and the precise question being asked. So what about making it an 
option in bipartite that you can either choose to divide by the realised links (give our 2007 paper as a ref, so people know it's comparable to that) or divide by the number of potential 
links, if that's the question people want to ask?''
There you go: it's your choice!
}

\value{
Depending on the selected indices, some or all of the below (returned as vector is \dQuote{degree distribution} was not requested, otherwise as list):
  \item{number of higher trophic species}{Number of species in the higher trophic level, i.e. \code{ncol(web)}.}
  \item{number of lower trophic species }{Number of species in the lower trophic level, i.e. \code{nrow(web)}.}
  \item{connectance}{Realised proportion of possible links (Dunne et al. 2002): sum of links divided by number of cells in the matrix (= number of higher times number of lower trophic level 
species). This is the standardised number of species combinations often used in co-occurrence analyses (Gotelli & Graves 1996)}  
  \item{web asymmetry}{Balance between numbers in the two levels: positive numbers indicate more lower-trophic level species, negative more higher-trophic level species; implemented as 
\code{(nrow(web)-ncol(web))/sum(dim(web))}; web asymmetry is a null-model for what one might expect in dependence asymmetry: see Bl¸thgen et al. (2007).}
  \item{links per species}{Mean number of links per species (qualitative): sum of links divided by number of species.}
  \item{number of compartments}{Compartments are sub-sets of the web which are not connected (through either higher or lower trophic level) to another compartment. Mathematically, 
they are Jordan blocks, but this implementation is rule-based (and fast). They are also nicely visualised in the \code{\link{visweb}} function.}
  \item{compartment diversity C.D.}{Shannon's diversity index across compartment sizes (i.e. number of participants); see Tylianakis et al. (2007).}
  \item{cluster coefficient}{The CC for a network is the average CC of its members. CC for each node, in turn, is simply the number of realised links devided by the number of possible links.
	 Introduced by Watts & Strogatz (1998) and described in Wikipedia under \url{http://en.wikipedia.org/w/index.php?title=Clustering_coefficient}. If you want to use Tore 
Opsahl's adaptation to two-modes, please see his function \code{clustering_tm} in \pkg{tnet}. To my knowledge, so far every study has used the ``wrong'' one, i.e.~the one used also here.}
  \item{weighted clustering coefficient}{When asking for ``weigthed cluster coefficient'', the (arithmetic mean) weighted version will automatically be returned unless the data are binary. 
The computation is based on \code{clustering_tm} in \pkg{tnet}. See there (and more on Tore Opsahl's webpages) for help.}
  \item{degree distribution}{Coefficients and fits for three different functions to degree distributions: exponential, power law and truncated power law. See \code{\link{degreedistr}} for details 
and references.}
  \item{mean number of shared hosts}{The simplest measure of co-occurrence and hence similarity in host preferences of the higher trophic level; based on Roberts & Stone (1990) and 
Stone & Roberts (1992).}
  \item{mean number of shared predators}{Same as before, only for the lower trophic level.}
  \item{togetherness}{Mean number of co-occupancies across all species-host-combinations; the whole matrix is scanned for submatrices of the form c(0,0,1,1), representing perfect 
matches of co-presences and co-absences. These are counted for each pairwise species combination, and averaged. Based on Stone & Roberts (1992).}
  \item{C score}{Mean (normalised) number of checkerboard combinations across all higher trophic level species. Values close to 1 indicate that there is evidence for disaggregation, 
e.g.~through competition. Value close to 0 indicate aggregation of species (i.e.~no repelling forces between species). See Stone and Roberts (1990) for details.}
  \item{V ratio}{Variance-ratio of species numbers to individual numbers within species for the higher trophic level. Values larger than 1 indicate positive aggregation, values between 
0 and 1 indicate disaggregation of species. See Schluter (1984) for details.}
  \item{discrepancy}{Discrepancy as proposed by Brualdi & Sanderson (1999); see \code{\link{discrepancy}} for details.}
  \item{nestedness}{Nestedness temperature of the matrix (0 means cold, i.e. high nestedness, 100 means hot, i.e. chaos). For details see \code{\link{nestedness}} and Rodriguez-Girones 
& Santamaría (2002). Notice that the function \code{\link{nestedness}} does not calculate any null model, simply because it is too computer-intensive. \code{networklevel} calls 
\code{nestedtemp}! If you are interested in the different null models, please use the function \code{\link{nestedness}} or \code{nestedtemp} in \pkg{vegan} directly.}
  \item{weighted nestedness}{A nestedness version that considers interaction frequencies (and is hence weighted), proposed by Galeano et al. (2007) and implemented in \code{\link{wine}}. 
It ranges between 1 (perfect nestedness) and 0 (perfect chaos). Note that this is the OPPOSITE interpretation of nestedness temperature!}
  \item{weighted NODF}{Another quantitative (=weighted) index for nestedness, building on NODF (see \code{nestednodf} in \pkg{vegan}). High values indicate nestedness. According to the analysis of Almeida-Neto et al. (2008, 2010), NODF is more consistent and ``better'' than usual measures of nestedness.}
  \item{dependence asymmetry/interaction strength asymmetry}{(selected using \option{index="ISA"}) Explaining dependence asymmetry will require more space than just a few lines. 
In essence, it is also a measure of specialisation, across both trophic levels. Proposed by Bascompte et al. (2006) and critised and alterations proposed by Blüthgen et al. (2007). 
The latter also show that dependence asymmetry can be almost entirely explained by web asymmetry (see above). Positive values (only possible of \option{DAmethod="Bluethgen"}) 
indicate higher dependence in the higher trophic level.}
  \item{specialisation asymmetry}{(selected using \option{index="SA"}) Asymmetry (higher vs. lower trophic level) of specialisation now based on d' (see \code{\link{dfun}}), which is
 insensitive to the dimensions of the web. Again, two options of calculation are available: the one proposed by Blüthgen et al. (2007), where they weight the specialisation value for each 
species by its abundance (\option{SAmethod="Bluethgen"}) or where d'-values are log-transformed (arguing that d'-values are indeed log-normally distributed: \option{SAmethod="log"}). 
Since the mean d-value for the lower trophic level is subtracted from that of the higher, positive values indicate a higher specialisation of the higher trophic level.}
  \item{extinction slope higher trophic level}{Slope of the secondary extinction sequence in the higher trophic level, following extermination of species in the lower trophic level; see details.}
  \item{extinction slope lower trophic level}{Slope of the secondary extinction sequence in the lower trophic level, following extermination of species in the higher trophic level; see details.}
  \item{robustness}{Calculates the area below the extinction curve; see \code{\link{robustness}} for details.}
  \item{niche overlap}{Mean similarity in interaction pattern between species of the same trophic level, calculated by default as Horn's index. Values near 0 indicate no common use of niches, 1 indicates perfect niche overlap. (In how far it makes sense for hosts of predators to call their commonness in enemies \dQuote{niche overlap} is a different issue. There are people calling predators \dQuote{negative resources} (couldn't be asked to look up the reference). I would probably rather call it similarity in predators.)}
  \item{generality}{Effective mean number of prey per predator; see Tylianakis et al. (2007).}
  \item{vulnerability}{Effective mean number of predator per prey; see Tylianakis et al. (2007).}
  \item{linkage density}{Mean number of interactions per species (quantitative); see Tylianakis et al. (2007).}
  \item{Fisher's alpha}{An alternative measure of interaction diversity (using \code{fisherfit} from \pkg{vegan}).}
  \item{HTL mean interaction diversity}{Mean diversity of prey per predator. HTL refers to Higher Trophic Level. HTLmid is analog to Shannon's diversity index for interactions (see below), 
but calculated for each species in the HTL, then averaged.}
  \item{LTL mean interaction diversity}{Mean diversity of predators per prey. LTL refers to Lower Trophic Level. LTLmid is analog to Shannon's diversity index for interactions (see below), 
but calculated for each species in the LTL, then averaged.}
  \item{interaction evenness}{Shannon's evenness for the web entries, treating zeros as no data. Note that the two options are rather different. Here a brief attempt to justify our default 
\code{"prod"}: By definition, IE = H/Hmax; H = -sum(p.i.mat*log(p.i.mat)), where p.i.mat = matrix/sum(entries in matrix). This means, when calculating H, we treat all possible links as 
species, and the interactions (cell values) as measure of their abundance. By definition, Hmax = ln(N). The key question is: What is the right value for N? Since we treat the matrix cells as 
species, it is (clearly?) the number of matrix cells, i.e. number of higher trophic level species x number of lower trophic level species. What else?
 Were we to use the interpretation of Tylianakis et al. (2007), then Hmax = ln(sum of entries in matrix). This means, we equate ``number of interactions'' (not links!) = ``number of species''. 
That means, each interaction is a species. What should that mean? Why should that measure ``interactions evenness''? Why should we move from a view of ``cells are species'' when 
calculating H to a view of ``interactions are species'' when calculating N? To say the least, it doesn't seem consistent. See details for more on this.}
  \item{Alatalo interaction evenness}{A different measure for web entry evenness, as proposed by M¸ller et al. (1999). }
  \item{Shannon diversity}{Shannon's diversity of interactions (i.e. network entries).}
  \item{H2}{H2' is a network-level measure of specialisation. It ranges between 0 (no specialisation) and 1 (complete specialisation). H2' is a measure of discrimination, i.e. calculated 
in comparison of no specialisation (see \code{\link{H2fun}} for details. To avoid confusion of keys (apostrophe vs. accent), we call the H2' only H2 here.}
}

\note{
All error and warning messages are (or at least should be) suppressed! If your web returns and NA for some of the indices, this can be because the index cannot be computed. For 
example, if the web is full (i.e. no 0-cells), extinction slopes cannot be fitted (singularity of gradient). Check if you can expect the index to be computable! If it is, and networklevel 
doesn't do it, let me know.

Some indices require rather long computation times on large webs. Particularly slow (for large webs) are: \dQuote{extinction slope} and \dQuote{nestedness}. For ordinary-sized webs 
(tens by tens), \code{networklevel} is reasonably fast. An example of a large web is \code{kato1990}: try it, possibly using \code{system.time} as a wrapper, and prepare to do some reading on the side ...

 If you want to increase the speed by omitting some indices, here a rough guide: Ask only for the indices you are interested in! Otherwise, here is the sequence of most time-consuming 
indices: 1. The slowest function is related to extinction slopes and robustness. Excluding both makes the function faster. 2. Degree distributions are somewhat time consuming. 
3. Fisher's alpha is computed iteratively and hence time consuming. 3. Nestedness and weighted nestedness are not the fastest of routines. 4. Number (and diversity) of compartments 
calls a recursive and hence relatively slow algorithm. 5. H2 and specialisation asymmetry require an iterative, heuristic search algorithm. Finally, excluding discrepancy can also 
moderately decrease computation time.
}

\references{ 
Almeida-Neto M, Loyola RD, Ulrich W, Guimaraes P, Guimaraes Jr. PR. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}, 1227--1239

Almeida-Neto M, Ulrich W (2010) A straightforward computational approach for measuring nestedness using quantitative matrices. \emph{Environmental Modelling & Software}, in press.

Bascompte, J., Jordano, P. and Olesen, J. M. 2006 Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. (2002) Quantitative descriptors of food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen N. 2007 Specialization, constraints and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Burgos, E., H. Ceva, R.P.J. Perazzo, M. Devoto, D. Medan, M. Zimmermann, and A. Maria Delbue (2007) Why nestedness in mutualistic networks? \emph{Journal of Theoretical Biology} \bold{249}, 307--313

Corso G, de Araújo AIL, de Almeida AM (2008) A new nestedness estimator in community networks. \emph{arXiv} 0803.0007v1 [physics.bio-ph]

Dormann, C.F., Fründ, J., Blüthgen, N., and Gruber, B. (2009) Indices, graphs and null models: analysing bipartite ecological networks. \emph{The Open Ecology Journal} \bold{2}, 7--24.

Dunne, J. A., R. J. Williams, and N. D. Martinez. 2002 Food-web structure and network theory: the role of connectance and size. \emph{Proceedings of the National Academy of Science USA} \bold{99}, 12917--12922

Galeano J, Pastor JM, Iriondo JM (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Gotelli, N. J., and G. R. Graves. 1996 Null Models in Ecology. Smithsonian Institution Press, Washington D.C.

Memmott, J., Waser, N. M. and Price M. V. 2004 Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Roberts, A. and Stone, L. 1990 Island-sharing by archipelago species. \emph{Oecologia} \bold{83}, 560--567

Rodríguez-Girónes M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935

Schluter, D. (1984) A variance test for detecting species associations, with some example applications. \emph{Ecology} \bold{65}, 998-1005.

Stone, L. and Roberts, A. (1990) The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79.

Stone, L. and Roberts, A. (1992) Competitive exclusion, or species aggregation? An aid in deciding. \emph{Oecologia} \bold{91}, 419--424

Tylianakis, J. M., Tscharntke, T. and Lewis, O.T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Watts, D. J. and Strogatz, S. (1998) Collective dynamics of \sQuote{small-world} networks. \emph{Nature} \bold{393}, 440--442
}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de}}

\seealso{ Some functions are implemented separately: \code{\link{H2fun}}, \code{\link{second.extinct}}, \code{\link{degreedistr}}, \code{\link{C.score}} and \code{\link{V.ratio}}
}

\examples{
data(Safariland)
networklevel(Safariland)
networklevel(Safariland, index="ALLBUTDD") #excludes degree distribution fits
}

\keyword{ package }