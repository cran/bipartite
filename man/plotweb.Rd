\name{plotweb}
\alias{plotweb}

\title{Visualize a bipartite interaction matrix (e.g. a foodweb)}
\description{
A two dimensional matrix is plotted as a bipartite graph.
}
\usage{
plotweb(web, method = "cca", empty = TRUE, labsize = 1, ybig = 1,
    y_width = 0.1, spacing = 0.05, arrow="no", col.interaction="grey80",
    col.pred = "grey10", col.prey="grey10", lab.space=1,
    lablength = NULL, sequence=NULL)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{method}{ Default method is \option{cca}, which leads to as few crossings of interactions as possible. The other option is \option{normal}, which leaves order as given by the matrix.}
  \item{empty}{logical; should empty columns or empty rows be omitted from plotting; defaults to true}
  \item{labsize}{ factor for size of labels, default is 1 }
  \item{ybig}{ vertical distance between upper and lower boxes, default is 1}
  \item{y_width}{ width of upper and lower boxes, default is 0.1 }
  \item{spacing}{ horizonatal distance between upper and lower boxes, default is 0.05}
  \item{arrow}{ display type of connection between upper and lower boxes, options are \option{up}, \option{down}, \option{both} and \option{no}, default is \option{no}, which is a polygonal connection between boxes. }
  \item{col.interaction}{color of interaction, default is grey80. }
  \item{col.pred}{color of upper boxes, default is grey10.}
  \item{col.prey}{color of lower boxes, default is grey10.}
  \item{lab.space}{sometimes it is neccessary to add additional space for labels below und above of the boxes, so all labes are shown, default is 1.}
  \item{lablength}{number of characters of labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{sequence}{list of two with two names vectors: \code{seq.pred} and \code{seq.prey}, which specify the order in which species are plotted. Cannot be set for \option{method="cca"}. Defaults to \code{NULL}, where the sequence remains as given or is determined by the CCA internally.}
}
\value{
  Returns a window with a bipartite graph of a foodweb.
}
\references{
Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. 2007 Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205
 }
\author{ Bernd Gruber \email{bernd.gruber@ufz.de} }

\seealso{ For a different plot of foodwebs see  \code{\link{visweb}} }
\examples{
data(Safariland)
plotweb(Safariland)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }


