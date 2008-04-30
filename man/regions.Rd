\name{regions}
\alias{regions}
\alias{regions<-}

\title{Retrieve regions data from cghRegions object.}
\description{
  This function accesses the regions values
  of assay data stored in an object derived from the \code{\link{cghRegions-class}}.
}
\usage{
regions(object)
regions(object) <- value
}
\arguments{
  \item{object}{Object derived from class \code{cghRegions}}
  \item{value}{Matrix with rows representing features and columns samples.}
}
\value{
  \code{regions} returns a matrix of regions values;
}
\author{Sjoerd Vosse}

\seealso{\code{\link{cghRegions-class}}}

\keyword{manip}
