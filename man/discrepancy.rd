\encoding{latin1}
\name{discrepancy}
\alias{discrepancy}

\title{ Calculates discrepancy of a matrix }

\description{
Discrepancy is the number of mismatches between a packed (binary) matrix and the maximally packed matrix (with same row sums)
}

\usage{
discrepancy(mat)
}

\arguments{
  \item{mat}{A matrix (or something that can be transformed into a matrix when \code{as.matrix} is called within the function) of species (in columns) on islands (in rows). If quantitative data are given (e.g. in a quantitative pollination network), these are internally transformed into a binary matrix.}
}

\details{
 Discrepancy is a way to measure the nestedness of a matrix. In a comparative study, Ulrich \& Gotelli (2007) showed discrepancy to outperform all other measures and hence recommend its use (together with a fixed-columns, fixed-rows null model, such as implemented in \code{commsimulator} in \pkg{vegan}, see example).

 This function follows the logic laid out by Brualdi \& Sanderson (1999), although, admittedly, I find their mathematical description highly confusing. Another implementation is given by the function \code{nesteddisc} in \pkg{vegan}. The reason to write a new function is simple: \code{nesteddisc} returns the wrong results, even when used on the simple example matrix given in Brualdi \& Sanderson. This is very strange, because Jari Oksanen is a very meticulous and efficient programmer, and you (the user) may want to trust his version more than this one. Having said that, values don't differ much between the two implementations. (The real reason why I wrote \code{discrepancy} anew was that I wasn't aware of \code{nesteddisc}. I was sitting on a train and I wanted to use this measure later on, so I put it into a function consulting only the orignal paper. When looking for the swap algorithm to create null models, which I somehow knew to exist in \pkg{vegan}, I stumbled across \code{nesteddisc}. If you are interested in the swap algorithm and come across this help page, let me re-direct you to \code{oecosimu} in \pkg{vegan}.)
 
 So what does it do: The matrix is sorted by marginal totals, yielding a matrix \bold{A}. Then, all 1s in \bold{A} are \dQuote{pushed} to the left to maximally compact the matrix, yielding \bold{P}. Discrepancy is now simply the number of disagreements between \bold{A} and \bold{P}, divided by two (to correct for the fact that every \dQuote{wrong} 1 will necessarily generate a \dQuote{wrong} 0).
}

\value{
Returns the number of mismatches, i.e. the discrepancy of the matrix from perfecct nestedness.
}

\references{ 
Brualdi, R.A. and Sanderson, J.G. (1999) Nested species subsets, gaps, and discrepancy. \emph{Oecologia} \bold{119}, 256--264

Ulrich, W. and Gotelli, N.J. (2007) Disentangling community patterns of nestedness and species co-occurrence. \emph{Oikos} \bold{116}, 2053--2061
}

\note{ Discrepancy is well-defined only for matrices that can be sorted uniquely. For matrices with ties no way to handle them has been proposed. For small matrices, or large matrices with many ties, this will lead to different discrepancy values. See also how \code{nesteddisc} in \pkg{vegan} handles this issue! (Thanks to Jari Oksanen for pointing this out!)
}

\author{ Carsten F. Dormann }

\seealso{\code{\link{nestedness}} for the most commonly used method to calculate nestedness, \code{\link{nestedness.corso}} for a new, unevaluated but very fast way to calculate nestedness; \code{nestedtemp} (another implementation of the same method used in our \code{nestedness}) and \code{nestedn0} (calculating the number of missing species, which has been shown to be a poor measure of nestedness) in \pkg{vegan}
}

\examples{
data(Safariland)
require(vegan)
nulls <- replicate(1000, discrepancy(commsimulator(Safariland, method="quasiswap")))
hist(nulls)
obs <- discrepancy(Safariland)
abline(v=obs, lwd=3, col="grey")
c("p value"=min(sum(nulls>obs), sum(nulls<obs))/length(nulls))
# calculate Brualdi & Sanderson's Na-value (i.e. the z-score):
c("N_a"=(unname(obs)-mean(nulls))/sd(nulls))
}

\keyword{package}
