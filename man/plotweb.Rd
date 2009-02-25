\encoding{latin1}
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
    bor.col.interaction ="black", bor.col.pred="black", bor.col.prey="black",
    lablength = NULL, sequence=NULL, low.abun=NULL, low.abun.col="green",
    bor.low.abun.col ="black", 
    high.abun=NULL, high.abun.col="red", bor.high.abun.col="black",
    text.rot=0)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{method}{ Default method is \option{cca}, which leads to as few crossings of interactions as possible. The other option is \option{normal}, which leaves order as given by the matrix.}
  \item{empty}{logical; should empty columns or empty rows be omitted from plotting; defaults to true}
  \item{labsize}{ factor for size of labels, default is 1 }
  \item{ybig}{ vertical distance between upper and lower boxes, default is 1}
  \item{y_width}{ width of upper and lower boxes, default is 0.1 }
  \item{spacing}{ horizonatal distance between upper and lower boxes, default is 0.05}
  \item{arrow}{ display type of connection between upper and lower boxes, options are \option{up}, \option{down}, \option{both}, \option{up.center}, \option{down.center}, \option{both.center} and \option{no}, default is \option{no}, which is a polygonal connection between boxes. }
  \item{col.interaction}{color of interaction, default is grey80. }
  \item{col.pred}{color of upper boxes, default is grey10.}
  \item{col.prey}{color of lower boxes, default is grey10.}
  \item{lab.space}{sometimes it is neccessary to add additional space for labels below und above of the boxes, so all labes are shown, default is 1.}
  \item{bor.col.interaction}{border color of interaction, default is black}
  \item{bor.col.pred}{border color of upper boxes, default is black}
  \item{bor.col.prey}{border color of lower boxes, default is black}
  \item{lablength}{number of characters of labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{sequence}{list of two with two names vectors: \code{seq.pred} and \code{seq.prey}, which specify the order in which species are plotted. Cannot be set for \option{method="cca"}. Defaults to \code{NULL}, where the sequence remains as given or is determined by the CCA internally.}
  \item{low.abun}{Vector with independent abundance estimates for the lower trophic level, NULL if none exists.}
  \item{low.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to green.}
  \item{bor.low.abun.col}{border color for depicting the abundance estimates for the lower trophic level, default is black}
  \item{high.abun}{Vector with independent abundance estimates for the higher trophic level, NULL if none exists.}
  \item{high.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to red.}
  \item{bor.high.abun.col}{border color for depicting the abundance estimates for the higher trophic level, default is black}
  \item{text.rot}{orientation of labels in the plot (to avoid overlapping of horizontal labels if dimension of the webs are high), default is 0 for horizontal labels, use text.rot=90 for vertical labels.}
}
\value{
  Returns a window with a bipartite graph of a foodweb.
}
\references{
Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205
 }
\author{ Bernd Gruber \email{bernd.gruber@ufz.de} }

\seealso{ For a different plot of foodwebs see  \code{\link{visweb}} }
\examples{
data(Safariland)
plotweb(Safariland)

# shorter labels
plotweb(Safariland, lablength=3, arrow="down")

# centered triangles for displaying interacions
plotweb(Safariland, text.rot=90, arrow="down.center", col.interaction="wheat2")

#orginal sequence, up arrows and different box width
plotweb(Safariland, method="normal", arrow="up", y_width=0.3, lablength=4)
# interactions as lines
plotweb(Safariland, arrow="both", y_width=0.05, text.rot=90, col.pred="blue", col.prey="green")

# add an abundance vector for lower trophic species 
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction="purple", y_width=0.05)

plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction ="red", bor.col.interaction="red", arrow="down")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }


