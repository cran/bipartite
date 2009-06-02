\encoding{latin1}

\name{swap.web}

\alias{swap.web}

\title{ Creates null model for bipartite networks}

\description{
Function to generate null model webs under the following constraints: 1. marginal 
totals are identical to those observed (as in \code{r2dtable}), 2. connectance is as observed (as in \code{shuffle.web}.)
}

\usage{swap.web(N, web, verbose=FALSE, c.crit=1e4)}

\arguments{
  \item{N}{ Number of desired nullmodel matrices. }
  \item{web}{An interaction matrix.}
  \item{verbose}{Should various verbal outputs of this function be shown? Defaults to FALSE, since it was
  mainly used during the debugging period.}
  \item{c.crit}{Sometimes the algorithm gets stuck in a very sparse matrix. Then c.crit sets the number of
  swaps it shall attempt before giving up and starting over on a new random matrix. Defauls to 10000.}
}

\details{
 This function is designed to behave similar to \code{r2dtable}, i.e. it returns a list of randomised matrices. 
 In addition to \code{r2dtable} is also keeps the connectance constant! 
 
 This function is thought of as a very constrained nullmodel for the analysis of bipartite webs. 
 It keeps two web properties constant: The marginal totals (as in \code{r2dtable} and the 
 number of links (and hence connectance). A comparison of \code{swap.web}- and \code{r2dtable}-based 
 webs allows to elucidate the effect of evolutionary specialisation, since the non-realised connections may
 represent \dQuote{forbidden links}. 
 
 This nullmodel is similar to the one employed by Vazquez et al. But while Vazquez starts by assigning 1s to the 
 allowed connections and then fills the web, \code{swap.web} starts with an \code{r2dtable}-web and successively 
 \dQuote{empties} it. The two approaches should result in very similar nullmodels, since both constrain marginal 
 totals and connectance.
 
 A few words on the way \code{swap.web} works. It starts with a random web created by \code{r2dtable}. Then, it finds,
 randomly, 2x2 submatrices with entries all larger than 0. Next, it subtracts the minimum value from the diagonal and
 adds it to the off-diagonal (minor diagonal). Thereby one cell becomes 0, but the column and row sums do not change.
 This idea is adapted from the swap-algorithm used in various binary null models by Nick Gotelli. If the random web has 
 too few 0s (which is I have yet to encounter), then the opposite strategy is applied.
 
 The algorithm in our implementation has some variations on finding the submatrix and constraining the number of 
 successless trials before starting on a new random matrix, but they are only for speeding up the process.
 
}

\value{
 A list of N randomised matrices with the same dimensions as the initial web.
}

\references{

 Vázquez, D. P., and M. A. Aizen (2003) Null model analyses of specialization in plant-pollinator interactions. \emph{Ecology} \bold{84} 2493-2501

 Vázquez, D. P., C. J. Melián, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin (2007) Species abundance and asymmetric interaction strength in ecological networks. \emph{Oikos} \bold{116}, 1120--1127

 For a very nice and thorough overview of nullmodels in general see:

 Gotelli, N. J., and G. R. Graves (1996) Null Models in Ecology. Smithsonian Institution Press, Washington D.C.
 

}

\author{ Carsten F. Dormann <carsten.dormann@ufz.de> }

\note{ 
  Long stories can be told about the swap algorithm. I am not the right person to do so, but for a much more detailed coverage of the topic, for many more ways to implement null models for \bold{binary} matrices, with various flavours of the swap and possible alternatives, first brew yourself a cup of tea and then check out the help pages of \code{commsimulator} in \pkg{vegan}. (As usual, Jari Oksanen has spend considerable care to implement even the most bizarr and abstruse way to move 0s and 1s around. His ecological advise between the lines make his package worthwhile already! I, personally, would use \code{method="quasiswap"}, as is done in the example to \code{\link{discrepancy}}.)

  \code{swap.web} is a very constraint nullmodel. You need to consider if it is the right one for your application!
}

\seealso{ \code{\link{r2dtable}} and \code{\link{shuffle.web}}}

\examples{

data(Safariland)

swap.web(Safariland, N=2)

}
\keyword{ package }

