\name{plot.cghRaw}
\alias{plot.cghRaw}
\alias{plot.cghSeg}
\alias{plot.cghCall}
\alias{plot.cghRegions}

\title{Plot aCGH data.}

\description{
  Please see the class descriptions for more details on the plot functions.
}

\usage{
  \method{plot}{cghRaw}(x, y, \dots)
  \method{plot}{cghSeg}(x, y, dotres=10, ylimit=c(-2,5), \dots)
  \method{plot}{cghCall}(x, y, dotres=10, ylimit=c(-5,5), \dots)
  \method{plot}{cghRegions}(x, y, \dots)
}

\arguments{
  \item{x}{ An object of class \code{\link{cghRaw}}, \code{\link{cghSeg}}, \code{\link{cghCall}}, or \code{\link{cghSeg}}.}
  \item{y}{ This argument is not used and should be missing.}
  \item{dotres}{ Every dotres-th log2-ratio is plotted. dotres=1 plots all data. However, higher values save a lot of space and 
     allow quicker browsing of the plots. }
  \item{ylimit}{ limits of the y-axis }
  \item{\dots}{ Arguments \code{plot}. }
}

\author{Sjoerd Vosse}

\seealso{\code{\link{cghRaw-class}}, \code{\link{cghSeg-class}}, \code{\link{cghCall-class}}, \code{\link{cghRegions-class}}}

