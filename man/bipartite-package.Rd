\encoding{latin1}

\name{bipartite-package}

\alias{bipartite-package}

\alias{bipartite}

\docType{package}

\title{
Analysis of bipartite ecological webs.
}

\description{
Bipartite provides functions to viualise webs and calculate a series of indices commonly used to describe pattern in ecological webs. It focusses on webs consisting of only two trophic levels, e.g. pollination webs or predator-prey-webs. Visualisation is important to get an idea of what we are actually looking at, while the indices summarise different aspects of the webs topology.}

\details{
We only had three types of bipartite webs in mind when writing this package: seed-disperser, plant-pollinator and host-parasitoid systems. In how far it makes sense to use these functionalities for other systems (or indeed for these systems) lies in the hands of the user. Please refer to the literature cited for details on the theory behind the indices.

Input for most analyses is an interaction matrix of m higher level species with n lower level species, i.e. an n x m matrix, where higher trophic level species are in columns, lower level in rows. Column and row names can be provided. This is fundamentally different from multi-dimensional webs, which are organised as k x k matrix, i.e. each species against each other. Such a format is incompatible with the functions we provide here.

The first step is to visualise the interaction web. Two functions are on offer here: one (\code{\link{visweb}}) simply plots the matrix in colours depicting the strength of an interaction and options for re-arranging columns and rows (e.g. to identify compartments or nesting). The other function (\code{\link{plotweb}}) plots the actual web with participants (as two rows of rectangles) connected by lines (proportional to interaction strength).

The second step is to calculate various indices describing network topography. There are two different levels this can be achieved at: the entire web (using function \code{\link{networklevel}}) or the individual species (using function \code{\link{specieslevel}}). All other functions in the package are helpers, although some can be called on their own and return the respective result (\code{\link{dfun}}, \code{\link{H2fun}} and \code{\link{second.extinct}} with \code{\link{slope.bipartite}}).

See function description for details and examples.

\tabular{ll}{
Package: \tab bipartite\cr
Type: \tab Package\cr
Version: \tab 1.02\cr
Date: \tab 2009-09-11\cr
License: \tab GPL \cr
}
}

\section{versionlog}{

\itemize{
  \item 1.02 (release date: 11-Sep-2009)
    \describe{
      \item{Function \code{\link{plotweb}}:}{Now more trophic networks can be plotted by staggering bipartite networks on top of each other. See multitrophic examples in \code{\link{plotweb}}.} 
      \item{Minor corrections to \code{\link{dfun}}:}{This function did not return ßemph{exactly} the values of the website-version. Jochen Fründ corrected this. Please read the help of the function (final paragraph) for details.}
      \item{New function \code{\link{nested}}:}{This convenience function collects the various ways to calculate nestedness of a network in order to facilitate comparison of nestedness analyses. To do so, it heavily borrows from \pkg{vegan}.}
      \item{Bug fix to \code{\link{discrepancy}}:}{Would return a silly value (half of the number of rows) for empty matrices. This had no effect when called by \code{\link{networklevel}}, since the matrix would have been emptied. Thanks to Roberto Molowny for reporting!}
      \item{Bug fix to \code{\link{networklevel}}:}{Due to a missing space, the option \option{ALLBUTDD} did not work properly. Thanks to Etienne Laliberté for reporting!}
    }
  \item 1.00 (release date: 06-Aug-2009)
    \describe{
      \item{Complete overhaul of \code{\link{networklevel}}:}{After a workshop on bipartite networks in ecology, a few more indices were added (Fisher's alpha diversity of interactions, mean interaction diversity, mean number of predators) and the whole output reorganised. It now follows a gradient from less to more interesting (in our view) indices, and from indices for binary to those for quantitative networks. Also, I added several options for which indices to report (\option{index="info", "binary", "quantitative", "topology"}). Most interestingly, perhaps, there is now a quantitative, Shannon-based series of indices. Starting with the \dQuote{mean interaction diversity} (i.e. the Shannon-diversity of interactions of a species, averaged across all species in that trophic level), over \dQuote{Shannon diversity} of interactions, to H2 (i.e. Shannon diversity scaled between max and min possible for this web characteristics).}
      \item{Additions to, and overhaul of \code{\link{specieslevel}}:}{Similar to the above item, some indices were added, and the output simplified when calling the option \option{index="ALLBUTD"} (only one D here!): a list with two matrices is now returned. Fisher's alpha for each species, vulnerability/generality and effective number of species for each target species are also now included. Index sequence has changed.}
      \item{New function \code{\link{PAC}}:}{Calculates the Potential for Apparent Competition following the formula in Müller et al. (1999) and Morris et al. (2005) and the suggestion by Becky Morris and Owen Lewis. More than a theoretical concept, it was experimentally shown to be relevant (Morris et al. 2004).}
      \item{Bug fix in \code{\link{H2fun}} and \code{\link{dfun}}:}{A line of code went missing at some point, so the maximum packing density was not optimal (but still good) in these functions. As a result, reported H2- and d-values were sometimes 0 when they should be only very close to 0.}
      \item{Switch for error reporting in \code{\link{degreedistr}}:}{By default now suppresses error reporting when \code{nls} fails to fit a degree distribution due to too few data points. This leaves the user of \code{\link{networklevel}} somewhat less confused. Warning message now also indicates for which trophic level there were too few data points.}
      \item{Bug fix in \code{\link{C.score}}:}{Failed when the web was very dense or very sparse, because the maximum number of checkerboard patterns was 0 then.}
    }
  \item 0.94 (release date: 01-Aug-2009)
    \describe{
      \item{New function \code{\link{nullmodel}}:}{A convenience wrapper function to generate different types of null models.}
      \item{Small changes to \code{\link{networklevel}}:}{This function returns a list of indices. If we exclude the computation of degree distribution fits, this would be coerced to a vector. We added the option \option{index="ALLBUTDD"} to calculate all indices BUT degree distributions. The output is then returned as vector. This is much more convenient when using networklevel on many data sets (using \code{sapply}).}
    }

  \item 0.93 (release date: 30-Jun-2009)
    \describe{
      \item{bug fix in \code{\link{slope.bipartite}} and \code{\link{robustness}}:}{The function always selected column 3 of the object, instead of 2 for lower and 3 for higher trophic level; thanks to Antonio Rivera for spotting and reporting this error! Notice that this error must have slipped in somewhere around version 0.90/0.91, because I checked and the results reported in the Open Ecology Journal paper are valid! I seem to remember that I (CFD) modified \code{\link{slope.bipartite}} when \code{\link{robustness}} was added: never change a wining team!}
    }
  \item 0.92 (release date: 02-Jun-2009)
    \describe{
      \item{more colour options in \code{\link{visweb}}:}{the arguments \option{box.border} and \option{box.col} now allow a specification of the colour of the boxes and their borders.}
      \item{bug fix \code{\link{empty}}:}{the function returned 0 for a 1x1 matrix. Although not written for such a case, it should still do what it says on the tin. Thanks to Mariano Devoto for spotting and reporting!}
      \item{bug fix \code{\link{wine}}:}{returned NA for square matrices.}
      \item{bug fix \code{\link{plot.wine}}:}{gave decimal places for row and column names for very small networks.}
    }

  \item 0.91 (release date: 06-May-2009)
  \describe{
    \item{new function \code{\link{wine}}:}{This function replaces the (now deprecated) function \code{nestedness.corso} in calculating a (weighted) nestedness for bipartite networks. It was developed and implemented by Marcelino de la Cruz, Juan M. Pastor, Javier Galeano and Jose Iriondo. It is also called by \code{\link{networklevel}}. A plotting function is also available, depicting the contributions of each observed link to the web's nestedness. - \code{nestedness.corso} is now removed from the package. The main reason is that it served as an interim solution for \code{\link{wine}}, and the Corso-way of calculating nestedness is just one more of already too many. If you intend to use \code{\link{wine}} on binary data and interpret that as a Corso-equivalent, beware of the following two main differences (thanks to Jose Iriondo for summarising them): \emph{First, the nestedness of Corso et al. varies between 0 and 1, with the highest nestedness is reached at 0 and 1 corresponds to random, whereas in \code{\link{wine}} is just the opposite (this is because the Manhattan distances are calculated with regard to opposite sides of the matrix). Secondly, the \sQuote{d} in the nestedness of Corso is the sums of the \sQuote{dij}s whereas in \code{\link{wine}}, the \sQuote{d} (= \sQuote{win} value of the object produced by \code{\link{wine}}) is the average of the \sQuote{dij}s above 0.} So, we recommend NOT to use wine for calculating Corso's nestedness, but rather download the source code for \code{nestedness.corso} from an older version of \pkg{bipartite}.}
    \item{new function \code{\link{robustness}}:}{A better way to quantify the effect of species loss on the extinctions in the other trophic level; kindly provided by Mariano Devoto. This index is also part of \code{\link{networklevel}}.}
    \item{new data set \code{ollerton2003}:}{Another quantitative pollination network from the NCEAS database (see \code{\link{ollerton2003}}).}
  }

  \item 0.90 (release date: 24-Mar-2009)
  \describe{
    \item{example \code{vazquez.example}:}{We added several new functionalities mainly to be able to analysis data and use network statistics as suggested by Vazquez et al. 2009. You can access this example by typing \code{?vazquez.example}.}
    \item{new feature in \code{\link{visweb}}:}{Can now plot different sized circles to represent interactions, as proposed by Vazquez et al. (2009).}
    \item{new feature in \code{\link{networklevel}}:}{We added an option to calculate interaction evenness either based on all possible links or just on realized links.}
    \item{new function \code{\link{sortweb}}:}{Can be used to sort webs in different ways.}
    \item{data set \code{inouye1988}:}{Another pollination network from the NCEAS database (see \code{\link{inouye1988}}).}
    \item{Function \code{\link{compart}}:}{We replaced the CA-based approach to detecting compartments by a comprehensive and recursive approach. The latter is not affected by ties in the data set (i.e. species with the same number of links). In quantitative webs and for the networks included in bipartite, the old function was working fine, but in more recent trials it failed to detect 2-species compartments. In turn, we had to adapt \code{\link{networklevel}} and \code{\link{plotweb}}.}
  }

	\item 0.85 (release date: 10-Mar-2009)
  \describe{
		\item{Function \code{\link{plotweb}}}{New feature: text labels can now be printed in different colours. All colours can be passed as vectors and vectors are recycled if not of appropriate length }
	}
	\item 0.84 (release date: 25-Feb-2009)
  \describe{
		\item{Function \code{\link{plotweb}}}{New feature: arrow=\dQuote{center.up}, arrow=\dQuote{center.down}, arrow=\dQuote{center.both}: this results in the standard presentation of bipartite networks, where interactions are displayed as centered triangles. }
	}

	\item 0.83 (release date: 9-Feb-2009)
	pdf of OEJ-paper added as vignette

	\item 0.82 (release date: 27-Jan-2009)
	\describe{
		\item{Function \code{\link{plotweb}}}{Two new features: 1. colors of borders of boxes and interactions can be specified via setting bor.col. 2. Labels can now be rotated by specifying text.rot=90.}
	}

	\item 0.81 (release date: 06-Jan-2009)
	\describe{
		\item{Function \code{\link{nodespec}} replaces the now obsolete function \code{functspec}.}{The inventor of functional specialisation, Bo Dalsgaard, understands the term \sQuote{functional specialisation} to be restricted to the way in which plants are specialised to pollinators. Although I disagree (and think that \sQuote{functional} actually means very little until defined for a given problem), I renamed the \sQuote{functional specialisation index} into \sQuote{node specialisation index}. This is not a particularly clever name, but at least it indicates that the position of nodes in a network is important when it is calculated.}
		\item{Bug in \code{\link{networklevel}}:}{Call to \code{nestedness.corso} overwrote results of \code{\link{nestedness}}.}
		\item{Various changes related to the analysis of very small webs.}{Very small webs should probably not be used at all for testing theories! Indices are usually VERY sensitive to the exact number of species, number of observations etc. Still, sometimes we simply want to calculate some index, and then \pkg{bipartite} should handle such small networks, too. So, when going through several dozens of very small webs (sometimes only containing one species in one of the two trophic levels), several functions did not perform correctly (usually to such minor programming issues such as matrices being converted to vectors when \code{[} was used with \option{drop=TRUE}, i.e. the default). Changes affected the functions \code{discrepancy}, \code{nestedness.corso}, \code{shuffle.web}, \code{compart} and \code{empty}. Their output remains identical, only they now also work for small webs.}
		\item{Changes to \code{nestedness.corso} with \option{weighted=TRUE}}{Galeano et al. do not describe how to deal with ties, nor do they make clear if the packed matrix should be sorted by marginal sum of links or marginal sum of interactions! Previously, we used marginal sum of interactions (because it is a weighted index), but now we moved to marginal sum of links, because that is how I interpret their paper after a third (or forth) re-reading. Furthermore, because the real maximum chaos cannot be derived (to my knowledge) algorithmically, we use the 95\% quantile of 500 randomisations as maximum. This will lead to a consistently overestimated nestedness, but it is less sensitive to the number of replicates than the max. (Also, there was an error in the description of the value returned: 0 is nested, 1 is chaos!)}
	}

	\item 0.8 (release date: 21-Dec-2008)
	\describe{
		\item{bug fix in \code{\link{shuffle.web}}}{Function didn't work correctly when the web contained more rows than columns. (Thanks to Anna Traveset for spotting and reporting!)}
		\item{minor modifications in networklevel}{Sometimes additional information is available and species can be included in a network, although they have no interactions with other species. In this case, one might want to use these species in the network, too. A new option (empty.web), allows to keep empty rows and columns, although for some functions an emptied web must be used (e.g. degree distributions).}
		\item{bug fix in \code{\link{H2fun}}}{As in the last \code{\link{H2fun}} bug fix, sometimes H2 became negative.}
		\item{bug (?) fix in \code{\link{networklevel}}}{Shannon diversity is based on the log of interactions. If this value is 0, as it is for most network entries, an NA is produced. As a consequence, Shannon's H (now also given in output) is based not on all interactions, but only those > 0. In consequence, Shannon evenness should also only be SH divided by the number of realised interactions (\code{log(sum(web>0))}). That is now the case.}

	}

	\item 0.74 (release date: 24-Oct-2008)
	\describe{
		\item{functional specialisation (\code{\link{functspec}}) bug fix}{Paths were double the true length, hence minimum was 2, rather than 1.}
		\item{\code{\link{H2fun}} bug fix}{Since the search for H2min is heuristic, H2uncorr can sometimes be lower than H2min; in that case, H2fun returned a value greater 1, while it should be one exactly.}
		\item{new function \code{nestedness.corso} }{Calculates (weighted) nestedness according to Corso et al. (2008) and Galeano et al. (2008).}
		\item{new function \code{\link{discrepancy}}}{Calculates discrepancy according to Brualdi & Sanderson (1999), deemed to be best ever measure of nestedness; also gives an example for a binary null model analysis based on \pkg{vegan}'s \code{commsimulator}.}
		\item{other}{Correction of several minor typos on the help pages; removal of "~" in help files; same citation style throughout; new cross references (especially for the nestedness functions); in \code{\link{networklevel}}, nestedness is now calculated using \pkg{vegan}'s \code{nestedtemp} due to matrix inversion problems reported for binmatnest.}
	}

	\item 0.73 (release date: 1-Sept-2008)
	\describe{	
		\item{new feature}{\code{\link{plotweb}}}{Named abundance-vector for each level can be used.}
		\item{new function \code{\link{plotweb2}} (not debugged!)}{For plotting tripartite networks.}
	}

	\item 0.72 (release date: 12-June-2008)
	\describe{
		\item{new function: functional specialisation with \code{\link{functspec}}}{See Dalgaard et al. (2008).}
		\item{new function: interface to \pkg{sna} through \code{\link{as.one.mode}}}{Allows calculation of path lengths, centrality, betweenness and other indices developed for one-mode networks.}
		\item{bug-fix}{Error in \code{\link{plotweb}} when no species labels were given.}
	}

}

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

Dormann, C.F., Fründ, J., Blüthgen, N., and Gruber, B. (2009) Indices, graphs and null models: analysing bipartite ecological networks. \emph{The Open Ecology Journal} \bold{2}, 7--24.

Galeano J., Pastor J.M. and Iriondo J.M. (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Memmott, J., Waser, N. M. and Price, M. V. (2004) Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2004. Experimental evidence for apparent competition in a tropical forest food web. \emph{Nature} \bold{428}, 310--313.

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2005. Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449--462.

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370.

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Vazquez, D. P. and Aizen, M. A. (2004) Asymmetric specialization: A pervasive feature of plant-pollinator interactions. \emph{Ecology} \bold{85}, 1251--1257

Vazquez, P.D., Chacoff, N.,P. and  Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology} in press.
}
\keyword{ package }

\examples{
data(Safariland)
plotweb(Safariland)
visweb(Safariland)
networklevel(Safariland)
specieslevel(Safariland)
}
