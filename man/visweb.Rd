\name{visweb}
\alias{visweb}

\title{ Plotfunction to visualize a bipartite foodweb }
\description{ this function draws a foodweb as a grid using a matrix. colnames and rownames are used as labels and entries in the matrix are visualized by text and colours

}

\usage{
visweb(web, type="nested",  prednames=TRUE, preynames=TRUE, labsize=1,
  plotsize=12, square="interaction", text="compartment", frame=NULL, textsize=1,
  textcol="red", pred.lablength=NULL, prey.lablength, clear=TRUE )
}

\arguments{
  \item{web}{ A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{type}{ type changes the sorting of rows and coloumns of the web and can be\cr
            \bold{nested  }: (sorted by row/colSums)\cr
            \bold{diagonal}: (highest number of interactions appear along the diagonal, good for showing compartments)\cr
            \bold{none    }: (as is)\cr
  }
  \item{prednames}{ labels can be switched of by prednames=F}
  \item{preynames}{ labels can be switched of by preynames=F}
  \item{labsize}{ factor for size of labels }
  \item{plotsize}{  size of plot (length of width or height), depending on the dimension of the web in cm,  default is 12 cm.}
  \item{square}{ square is used to indicate number of interactions  and belonging to compartments by coloured grid cells\cr
            \bold{interaction}: (level of grey indicates the number of intercation, white means no interaction)\cr
            \bold{compartment}: (level of grey indicates belonging to the same compartment)\cr
            \bold{black      }: (black grid cells if number of interaction is bigger than one)\cr
            \bold{none       }: (no coloured squares are drawn)\cr
  }
  \item{text}{ number of interactions or belonging are plotted into each grid cell\cr 
            \bold{interaction}: (number of interactions are drawn)\cr
            \bold{compartment}: (belonging to same compartment indicated by capitel Letters)\cr
            \bold{none       }: (no text is shown)\cr
  }
  \item{frame}{ a frame is drawn around each compartment (frame=T), best used with type=\dQuote{diagonal}  }
  \item{textsize}{ factor for size of text in squares, default=1 }
  \item{textcol}{  color of text in grid cells, default =\dQuote{red}}
  \item{pred.lablength}{length of predators (upper) labels that should be displayed}
  \item{prey.lablength}{length of prey (lower) labels that should be displayed}
  
  \item{clear}{ delete species with no interactions (compulsory done for \dQuote{nested} and \dQuote{diagonal}) }
}
\details{
}
\value{
  A plot window with appropriate size according to the dimensions of the web.
}
\references{ ... have to find appropriate ones  }
\author{ Bernd Gruber }

\seealso{ For a different plot on foodweb see \code{\link{plotweb}} }
\examples{

data(Safariland)
visweb(Safariland)
visweb(Safariland, type="diagonal", square="compartment", text="none", frame=TRUE)
}
\keyword{ package }

