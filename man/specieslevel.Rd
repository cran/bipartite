\encoding{latin1}
\name{specieslevel}

\alias{specieslevel}

\title{Calculate various indices for network properties at the species level}

\description{
Apart from the properties of the entire web, also its participants can be described
specifically. Various simple numbers and indices are calculated and returned.
}

\usage{
specieslevel(web, index="ALL", logbase="e", low.abun=NULL, high.abun=NULL)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher
  trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{index}{Vector of indices to be calculated for each trophic level of the web;
            options are: \dQuote{species number} for number of species, 
            \dQuote{degree},  \dQuote{dependence}, \dQuote{interaction} for interaction push/pull (our version of dependence asymmetry: see details), 
            \dQuote{PSI} for pollination service index (or pollinator
            support index, depending on the trophic level), \dQuote{NS} for node
            specialisation index, \dQuote{Fisher} for Fisher's alpha index, \dQuote{diversity} for Shannon diversity of interactions of that species, \dQuote{effective partners} for the effective number of interacting partners (e.g. predators/pollinators or hosts/plants),  \dQuote{d} for Blüthgen's d, \dQuote{strength} as sum of dependencies or \dQuote{ALL} for all the aforementioned (default). \dQuote{ALLBUTD} excludes only the dependence matrix and leads to the output being simplified to a matrix per trophic level.}
  \item{logbase}{Shall the various diversity indices be calculated to the base of e (default) or 2? Log2 is the proposal for generality and vulnerability by Bersier et al. (2002), while Shannon uses ln. The choice of the base will not affect the results qualitatively, only by a scaling factor.}
  \item{low.abun}{Optional vector of independent abundances of the \emph{lower} trophic level to be passed on to \code{\link{dfun}}.}
  \item{high.abun}{Optional vector of independent abundances of the \emph{higher} trophic level to be passed on to \code{\link{dfun}}.}
}

\details{
This function implements a variety of the many (and still procreating) indices describing
species properties. Some are embaracingly simple (such as number of interacting species for each species).
Others are variations on Shannon's diversity index applied to within species.

Most indices are straightforward, one-line formulae; some, such as d', also require a re-arranging of the matrix. We (Dormann, Blüthgen, Gruber) came up with a new one, called \dQuote{Pollination Service Index} or psi, for which a few more details seem appropriate.

%PSI calculates the average contribution per visit for each pollinator species (which in itself depends
%on the specialisation and abundance of the bees, as well as the abundance of the plant species). It is hence a measure
%of the \sQuote{usefullness} of each pollinator for the observed plant community. The value is ranged
%between 0 (pollinator largely useless for this plant community) and 1 (perfect specialisation by abundant pollinators).
%The usefulness of this value lies in the fact that a secondary extinctions will increase as the psi-value of a
%extinct pollinator increases: psi=1 implies a full dependence of at least one species on this pollinator, while a
%psi of 0 indicates a pollinators redundancy.

PSI is comprised of three calculation steps: firstly, we calculate, for each pollinator species, the proportion to
which it visits each plant species (or, phrased anthromorphically, the number to the question: To which proportion do I
visit dandelion?). Secondly, we calculate the proportion to which a plant is visited by each bee species (Which
proportion of my pollinators are red mason bees?). Multiplying, these two proportions gives the portion of own pollen
for each plant species (because this depends both on a pollinators specialisation (step 1) and the plant's
specific receptiveness (step 2). Finally, we sum the proportions own pollen delivered across all plant species.
This value is the PSI-value. At its maximum, 1, it shows that all pollen delivered to exactly one plant species is
derived from one monolectic pollinator. At its minimum, 0, it indicates that a pollinator is irrelevant to all plant
species. Note that PSI is independent of the abundance of plants or pollinators: a bee been found only once on a
plant species visited by no-one else receives a PSI of 1, even if in total 14 million visits were recorded.

We envisage a penalty for the fact that a pollinator has to make two (more or less successive) visits to the same
plant species: the first to take the pollen up, the second to pollinate the next. Thus, using 2 as an exponent
in step 1 would simulate that a pollinator deposits all pollen at every second visit (or flies to their nest). That is
certainly a too strong penalisation. At present we set the exponent to 1, because the step of controlling for
\dQuote{pollen purity} is already a major improvement. Also, we have no idea to which extent pollen gets mixed
and/or lost during foraging flights, and the true exponent remains elusive.

}
\value{
For both the \dQuote{higher trophic level} and the \dQuote{higher trophic level} a list with the following components:

  \item{number of species}{Sic.}
  \item{species degree}{Sum of interactions per species.}
  \item{dependence}{Web entry divided by column totals (for higher trophic level) or row totals (for lower
        trophic level). Expresses the proportion of interactions with each species, and sums to 1 for each.}
  \item{strength}{Sum of dependencies of each species (proposed in Bascompte et al. 2006).}
  \item{interaction push/pull}{Direction of interaction asymmetry: positive values indicate that a species affects the species of the other level it interacts with stronger then they affect it (\dQuote{pusher}); negative values indicate that a species is, on average, on the receiving end of the stick (\dQuote{being pulled}); formula based on Vazquez et al (2007), but push/pull is our own nomenclature. Values are highly correlated with	species strengths (see below), but standardised to fall between -1 (being pulled) and 1 (pushing).}
  \item{Pollination Service Index PSI}{for the higher trophic level, and the equivalent Pollinator Support Index for the lower trophic level. Expresses the value of, say, a pollinator for all plant species as the sum (across all plant species) of the proportion of pollen deposited per visit. The more specialised the pollinator and the fewer pollinators the plants it pollinates, the higher is its value for the plant community. It is a really cool measure, but it's not easy to get one's head around it. Alternatively, one can see it as the inverse of redundancy of the species for the web: the higher the value, the less redundant.}
  \item{node specialisation index}{Another measure of specialisation, based on the path length between any two higher-trophic level species. Species sharing hosts/prey have an FS-value of 1. See specific function \code{\link{nodespec}} for details, problems and reference.}
  \item{Fisher alpha}{Fisher's alpha diversity for each species (see \code{fisher.alpha} in \pkg{vegan} for details).}
  \item{Partner diversity}{Shannon diversity (when using \option{logbase="e"}) or per-species generality/vulnerability (when using \option{logbase=2}) of the interactions of each species. See also \code{\link{networklevel}} for the aggregated version of this index (i.e. averaged across all species in a trophic level).}
  \item{effective partners}{\option{logbase} to the power of \dQuote{partner.diversity}: Bersier et al. (2002) interpret this as the effective number of partners, if each partner was equally common. Note that \dQuote{partner} is a bit euphemistic when it comes to predator-prey or host-parasitoid networks.}
  \item{d}{Specialisation of each species based on its discrimination from random selection of
        partners. More specifically, it returns d', which is calculated based on the
        raw \code{d}, \code{dmin} and \code{dmax} for each species (see \code{\link{dfun}}. See Blüthgen et al. (2006) for details.}
}

\references{
  Bascompte, J., Jordano, P. and Olesen, J. M. (2006) Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

  Blüthgen, N., Menzel, F. and Blüthgen, N. (2006) Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 9

  Vázquez, D. P., Melián, C. J., Williams, N. M., Blüthgen N., Krasnov B. R. and Poulin, R. (2007) Species abundance and asymmetric interaction strength in ecological networks. \emph{Oikos} \bold{116}, 1120--1127
}

\author{ Carsten F. Dormann \email{carsten.dormann@ufz.de}}

\note{
Dependencies are still an open field of debate. Dependencies are calculated as the value in a matrix divided by the \code{rowSums} (for the lower trophic level) or the \code{colSums} (for the higher trophic level). As such, any pollinator observed only once will receive a dependency-value of 1, indicating perfect dependence on this plant species. That may or may not be true. In any case it is based on a sample size of 1, that is why the dependency asymmetry (which is based on the dependencies for both trophic levels) has come under criticism and may be rather sensitive to singletons.

We here provide the code to calculate the strength of a species (i.e. sum of its dependencies), based on the current proposal by Bascompte et al. (2006). It may be a good idea to remove all singletons from the web before calculating this index, to investigate whether it is indeed driven by those scarce observations.

The maximum value of the uncorrected d is not a trivial issue. We treat it here in the same way as given in the BMC Ecology paper, but please have a look at the (raw) code for further comments.
}

\seealso{ \code{\link{networklevel}} for some further comments; \code{\link{dfun}}, \code{\link{nodespec}}, which are called by this function }

\examples{
data(Safariland)
specieslevel(Safariland)
specieslevel(Safariland, index="ALLBUTD")[[2]]
}


\keyword{ package}
