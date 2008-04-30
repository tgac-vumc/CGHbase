\name{avedist}
\alias{avedist}
\alias{nclone}

\title{Retrieve regions information from cghRegions object.}
\description{
  This function accesses the regions information
  stored in the featureData of an object derived from the \code{\link{cghRegions-class}}.
}
\usage{
avedist(object)
nclone(object)
}
\arguments{
  \item{object}{Object derived from class \code{cghRegions}}
}
\value{
  \code{avedist} returns a vector containing the Average L1-distance of clone signatures to the medoid signature;
  \code{nclone} returns a vector containing the number of clones that is included in each region;
}
\author{Sjoerd Vosse}

\seealso{\code{\link{cghRegions-class}}}

\keyword{manip}
