\name{frequencyPlot}
\alias{frequencyPlot}
\title{ Visualization of aCGH regions. }
\description{
  This function creates a frequency plot for aCGH regions.
}
\usage{
	frequencyPlot(x, y, ...)
}
\arguments{
  \item{x}{ An object of class \code{\link{cghRegions}}.}
  \item{y}{ This argument is not used and should be missing. }
  \item{\dots}{ Arguments \code{plot}. }
}
\details{
  We find plotted on the x-axis the array probes sorted by chromosomal position. The vertical bars represent the frequency of gains and losses across your samples. The black bars represent gains, the gray bars represent losses.  
}
\value{
  This function creates a plot.
}
\references{
	Mark A. van de Wiel and Wessel N. van Wieringen (2007). CGHregions: Dimension Reduction for Array CGH Data
	with Minimal Information Loss. \emph{Cancer Informatics, 2}, 55-63.
}
\author{ Mark van de Wiel and Sjoerd Vosse}
\examples{
	\dontrun{
		data(WiltingRegions)
		frequencyPlot(WiltingRegions)
	}
}
\keyword{ misc }
