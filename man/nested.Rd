\encoding{UTF-8}

\name{nested}
\alias{nested}

\title{
Calculates any of several measures of nestedness 
}

\description{
Wrapper function calling one, several or all currently implemented nestedness measures
}

\usage{
nested(web, method = "binmatnest2", ..., rescale = FALSE)
}

\arguments{
  \item{web}{A matrix with elements of a set (e.g., plants) as rows, elements of a second set (e.g., pollinators) as columns and number of interactions as entries.}
  \item{method}{One or more of the following: discrepancy, discrepancy2, binmatnest, binmatnest2, NODF, NODF2, C.score, checker, , nested NODF, wine, ALL. See details for details on each method}
  \item{\dots}{Arguments passed on to other nestedness functions. Options need to be specified (i.e. no positional parsing).}
  \item{rescale}{Should all measures be rescaled so that higher values mean higher nestedness? Defaults to FALSE, i.e. the standard interpretation of each measure is maintained.}
}

\details{
  There are five different measures currently available:
  \describe{
    \item{1}{binmatnest calculates nestedness temperature following the function \code{\link{nestedness}} (0 = cold = highly nested; 100 = hot = not nested at all). It uses the original 
program of Miguel Rodriguez-Girones, only called from R; binmatnest2, in contrast, is the implementation in \code{nestedtemp} of the same algorithm by Jari Oksanen. Because 
binmatnest sometimes (and to us unexplicably) invert the matrix, we prefer the binmatnest2 option.}
    \item{2}{Discrepancy calculates the number of non-nested 0s and 1s. While \code{discrepancy} calls the function with the same name, discrepancy2 calls \code{nesteddisc}, 
which handles ties differently. Most of the time, these two should deliver very, very similar results. Higher values indicate lower nestedness.}
    \item{3}{NODF is the nestedness measure proposed by Almeida-Neto et al., correcting for matrix fill and matrix dimensions. Values of 0 indicate non-nestedness, those of 100 perfect nesting. NODF2 sorts the matrix before calculating the measure. NODF is, I understand, closer to the version presented in the paper, while NODF2 seems to make more sense for comparisons across different networks (because it is independent of the initial presentation of the matrix). Both call \code{nestednodf} in \pkg{vegan}. (Yes,  I initially programmed NODF myself, only to find that it was there already. Luckily, there was a perfect agreement between my (depricated) version and nestednodf.) A weighted version is also now available, following the paper of Almeida-Neto and Ulrich (2010).}
    \item{4}{\code{\link{C.score}} calculates the number of checkerboard pattern in the matrix. As default, it normalises this value between min and max, so that values of 0 
indicate no checkerboards (i.e. nesting), while a value of 1 indicates a perfect checkerboard. checker is the non-normalised version, based on \code{nestedchecker}.}
    \item{5}{wine is the only nestedness measure using the information on the weight of a link. See \code{\link{wine}} for details.}
  }


}
\value{
A vector with values for each of the selected nestedness measures.
}

\references{
Almeida-Neto, M., Guimaraes, P., Guimaraes, P.R., Loyola, R.D. and Ulrich, W. 2008. A consistent metric for nestedness analysis in
     ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}, 1227--1239.

Blüthgen, N., J. Fründ, D. P. Vazquez, and F. Menzel. 2008. What do interaction network metrics tell us about specialisation and biological traits? \emph{Ecology} \bold{89}, 3387--3399.

Brualdi, R.A. and Sanderson, J.G. 1999. Nested species subsets, gaps, and discrepancy. \emph{Oecologia}  \bold{119}, 256--264.
     
Galeano, J., Pastor, J.M., Iriondo and J.M. 2008. Weighted-Interaction  Nestedness Estimator (WINE): A new  estimator to calculate over
     frequency matrices. \emph{arXiv} 0808.3397v2 [physics.bio-ph]
     
Rodríguez-Gironés, M.A.  and Santamaría, L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence
     matrices. \emph{J. Biogeogr.}  \bold{33}, 924--935.

Stone, L. and Roberts, A. 1990. The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79.

Almeida-Neto, M. and Ulrich, W. 2010. A straightforward computational approach for measuring nestedness using quantitative matrices. Environmental Modelling & Software, in press.
}

\author{
Carsten F. Dormann \email{carsten.dormann@ufz.de}
}

\note{
The idea behind this function is to encourage the comparison of different nestedness measures. That does not mean, we necessarily see much ecological sense in them (see, e.g., the paper by Bl¸thgen et al. 2008).
}


\seealso{
\code{\link{C.score}}, \code{\link{wine}}, \code{\link{nestedness}}, \code{\link{discrepancy}}; and, within \pkg{vegan}: \code{nestedtemp},   \code{nestedchecker}, \code{nesteddisc}, \code{nestednodf}
}
\examples{
\dontrun{
data(Safariland)
nested(Safariland, "ALL")
nested(Safariland, "ALL", rescale=TRUE)
# illustration that non-normalised C.score and checker are the same:
nested(Safariland, c("C.score", "checker"), normalise=FALSE)
}
}

\keyword{ package }

