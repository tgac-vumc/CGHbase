\name{cghCall}
\docType{class}
\alias{class:cghCall}
\alias{cghCall}

\alias{cghCall-class}

\alias{copynumber,cghCall-method}
\alias{copynumber<-,cghCall,matrix-method}
\alias{segmented,cghCall-method}
\alias{segmented<-,cghCall,matrix-method}
\alias{calls,cghCall-method}
\alias{calls<-,cghCall,matrix-method}
\alias{probloss,cghCall-method}
\alias{probloss<-,cghCall,matrix-method}
\alias{probnorm,cghCall-method}
\alias{probnorm<-,cghCall,matrix-method}
\alias{probgain,cghCall-method}
\alias{probgain<-,cghCall,matrix-method}
\alias{probamp,cghCall-method}
\alias{probamp<-,cghCall,matrix-method}

\alias{chromosomes,cghCall-method}
\alias{bpstart,cghCall-method}
\alias{bpend,cghCall-method}

\alias{initialize,cghCall-method}
\alias{plot.cghCall,cghCall,missing-method}
\alias{summaryPlot,cghCall,missing-method}

\title{
  Class to contain and describe called array comparative genomic hybridization data.
}

\description{
  Container for aCGH data and experimental
  metadata. \code{cghCall} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires the following matrices of equal dimension
  as assayData members:
  \itemize{
	\item \code{copynumber} Log2 copynumber ratios.
	\item \code{segmented}	Segmented log2 ratios.
	\item \code{calls}		  Called copynumer values.
	\item \code{probloss}		Loss probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probnorm}		Normal probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probgain}		Gain probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probamp}		Optional amplification probabilities as returned by \code{\link[CGHcall]{CGHcall}} when run with \code{nclass=4}.
  }
  Furthermore, columns named \code{Chromosome}, \code{Start}, and \code{End} are 
  required as featureData members, containing feature position information.
}
\section{Extends}{
   Directly extends class \code{\link[Biobase:class.eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('cghCall',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    copynumber = [matrix],
    segmented = [matrix],
    calls = [matrix],
    probloss = [matrix],
    probnorm = [matrix],
    probgain = [matrix],
    featureData = [AnnotatedDataFrame],
    ...)
  }

  An object of class \code{cghCall} is generally obtained as output
  from \code{\link[CGHcall]{CGHcall}}.

}
\section{Slots}{
  Inherited from \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices with equal
	dimensions, and with column number equal to
	\code{nrow(phenoData)}. \code{assayData} must contain the following matrices
  \itemize{
	\item \code{copynumber} Log2 copynumber ratios.
	\item \code{segmented}	Segmented log2 ratios.
	\item \code{calls}		  Called copynumer values.
	\item \code{probloss}		Loss probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probnorm}		Normal probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probgain}		Gain probabilities as returned by \code{\link[CGHcall]{CGHcall}}.
	\item \code{probamp}		Optional amplification probabilities as returned by \code{\link[CGHcall]{CGHcall}} when run with \code{nclass=4}.
  }	
  with rows represening array probes
	and columns representing samples. Additional matrices of
	identical size (e.g., representing measurement errors) may
	also be included in \code{assayData}. Class:\code{\link[Biobase:class.AssayData]{AssayData-class}}}
      \item{\code{phenoData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{featureData}:}{An \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} with columns
	\code{Chromosome}, \code{Start}, and \code{End} containing array element position 
	data.}
      \item{\code{experimentData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{
     \item{\code{copynumber(cghCall)}, \code{copynumber(cghCall,matrix)<-}}{Access and
       set elements named \code{copynumber} in the \code{AssayData-class}
       slot.}
     \item{\code{segmented(cghCall)}, \code{segmented(cghCall,matrix)<-}}{Access and
       set elements named \code{segmented} in the \code{AssayData-class}
       slot.}
     \item{\code{calls(cghCall)}, \code{calls(cghCall,matrix)<-}}{Access and
       set elements named \code{calls} in the \code{AssayData-class}
       slot.}
     \item{\code{probloss(cghCall)}, \code{probloss(cghCall,matrix)<-}}{Access and
       set elements named \code{probloss} in the \code{AssayData-class}
       slot.}
     \item{\code{probnorm(cghCall)}, \code{probnorm(cghCall,matrix)<-}}{Access and
       set elements named \code{probnorm} in the \code{AssayData-class}
       slot.}
     \item{\code{probgain(cghCall)}, \code{probgain(cghCall,matrix)<-}}{Access and
       set elements named \code{probgain} in the \code{AssayData-class}
       slot.}
     \item{\code{chromosomes}, \code{bpstart}, \code{bpend}}{Access the chromosomal positions stored in \code{featureData}}
      \item{plot}{Create a plot containing log2ratios, segments and call probabilities ordered by chromosomal 
     position. TWO EXTRA OPTIONS PLUS DEFAULTS: 
     dotres=10. Every dotres-th log2-ratio is plotted. dotres=1 plots all data. However, higher values save a lot of space and 
     allow quicker browsing of the plots. ylimit=c(-5,5): limits of the y-axis}
     \item{summaryPlot}{Create a plot summarizing the call probabilities of all samples}
    }

   See \code{\link[Biobase:class.eSet]{eSet}} for derived methods.   
}

\author{Sjoerd Vosse}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet-class}}, \code{\link{cghRaw-class}}, \code{\link{cghSeg-class}}
}

\examples{

	# create an instance of cghCall
	new("cghCall")
	
	# load an instance of cghCall
	data(WiltingCalled)
  
  # plot the first sample
  plot(WiltingCalled[,1])
  
  # plot the first chromosome of the first sample
  plot(WiltingCalled[chromosomes(WiltingCalled)==1,1])

  # get the copynumber values of the third and fourth sample
  log2ratios <- copynumber(WiltingCalled[,3:4])

  # get the names of the samples
  sampleNames(WiltingCalled)

  # get the names of the array elements
  featureNames(WiltingCalled)

}
\keyword{classes}
