\name{cghRegions}
\docType{class}
\alias{class:cghRegions}
\alias{cghRegions}

\alias{cghRegions-class}

\alias{regions,cghRegions-method}
\alias{regions<-,cghRegions,matrix-method}

\alias{chromosomes,cghRegions-method}
\alias{bpstart,cghRegions-method}
\alias{bpend,cghRegions-method}
\alias{nclone,cghRegions-method}
\alias{avedist,cghRegions-method}

\alias{initialize,cghRegions-method}
\alias{plot.cghRegions,cghRegions,missing-method}
\alias{frequencyPlot,cghRegions,missing-method}

\title{
  Class to contain and describe array comparative genomic hybridization regions data.
}

\description{
  Container for aCGH regions data and experimental
  metadata. \code{cghRegions} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires a matrix named \code{regions} as
  assayData member. Furthermore, columns named \code{Chromosome}, \code{Start}, \code{End},
  \code{Nclone}, and \code{Avedist} are required as featureData members, containing region 
  and position information.
}
\section{Extends}{
   Directly extends class \code{\link[Biobase:class.eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('cghRegions',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    regions = [matrix],
    featureData = [AnnotatedDataFrame],
    ...)
  }

	An object of this class is generally obtained by running the function \code{\link[CGHregions]{CGHregions}}.
}
\section{Slots}{
  Inherited from \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices with equal
    dimensions, and with column number equal to
    \code{nrow(phenoData)}. \code{assayData} must contain a matrix
    \code{regions} with rows represening regions
    and columns representing samples. Additional matrices of
    identical size (e.g., representing measurement errors) may
    also be included in \code{assayData}. Class:\code{\link[Biobase:class.AssayData]{AssayData}}}
      \item{\code{phenoData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{featureData}:}{An \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} with columns
    \code{Chromosome}, \code{Start}, \code{End}, \code{Nclone}, and \code{Avedist} containing region 
    and position information.}
      \item{\code{experimentData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{
     \item{\code{regions(cghRegions)}, \code{regions(cghRegions,matrix)<-}}{Access and
       set elements named \code{regions} in the \code{AssayData-class}
       slot.}
     \item{\code{chromosomes}, \code{bpstart}, \code{bpend}, \code{nclone}, \code{avedist}}{Access the region and position
      information stored in \code{featureData}}
     \item{plot.cghRegions}{Create a plot displaying chromosomes on the Y-axis and base pair position on the X-axis. A new region is displayed by a slight jump with respect to the previous region. Each region is displayed as a bi-colored segment, the lower and upper part of which correspond to the proportions pl and pg of samples with a loss (red) or gain (green), respectively. The color coding is displayed as well: 1: pl (pg) < 10\%; 2: 10\% = pl (pg) < 30\%; 3:30\% = pl (pg) < 50\%; 4: pl (pg) = 50\%.}
     \item{frequencyPlot}{Create a frequency plot}
    }

   See \code{\link[Biobase:class.eSet]{eSet}} for derived methods. Annotation functionality is not yet supported.   
}

\author{Sjoerd Vosse}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet}}, \code{\link{cghRaw-class}}, \code{\link{cghSeg-class}}, \code{\link{cghCall-class}}
}

\examples{

	# create an instance of cghRegions
	new("cghRegions")
	
	# load an instance of cghRegions
	data(WiltingRegions)
		
	# plot all region data
	plot(WiltingRegions)
	# make a frequency plot
	frequencyPlot(WiltingRegions)
	
	# extract the region values
	values <- regions(WiltingRegions)
	
	# get the names of the samples
	sampleNames(WiltingRegions)

}
\keyword{classes}
