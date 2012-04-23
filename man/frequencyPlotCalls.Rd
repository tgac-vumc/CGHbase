\name{frequencyPlotCalls}
\alias{frequencyPlotCalls}
\title{ Visualization of aCGH profiles. }
\description{
  This function creates a frequency plot for aCGH profiles.
}
\usage{
	frequencyPlotCalls(x,main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37', ...)
}
\arguments{
  \item{x}{ An object of class \code{\link{cghCall}}.}
  \item{main}{Title of plot}
  \item{gaincol}{Color to use for gains}
  \item{losscol}{Color to use for losses}
  \item{misscol}{Missings}
  \item{build}{Build of Humane Genome.Either \code{GRCh37}, \code{GRCh36}, \code{GRCh35} or \code{GRCh34}}
  \item{\dots}{ Arguments \code{plot}. }
}
\details{
  We find plotted on the x-axis the array probes sorted by chromosomal position. 
  The vertical bars represent the frequency of gains or losses. 
}
\value{
  This function creates a plot.
}
\references{
  Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. \emph{Bioinformatics, 23}, 892-894.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
	 \dontrun{
   data(Wilting)
   rawcgh <- make_cghSeg(Wilting)
   normalized <- normalize(rawcgh)
   segmented <- segmentData(normalized)
   called <- CGHcall(segmented,cellularity= rep(0.75, 3))
   frequencyPlotCalls(called)
   }
}
\keyword{ misc }
