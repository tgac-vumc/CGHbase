### cghRaw

setGeneric("plot.cghRaw",   function(x, y, ...)     standardGeneric("plot.cghRaw"))
setGeneric("copynumber",    function(object)        standardGeneric("copynumber"))
setGeneric("copynumber<-",  function(object, value) standardGeneric("copynumber<-"))
setGeneric("chromosomes",   function(object)        standardGeneric("chromosomes"))
setGeneric("bpstart",       function(object)        standardGeneric("bpstart"))
setGeneric("bpend",         function(object)        standardGeneric("bpend"))

### cghSeg

setGeneric("plot.cghSeg",   function(x, y, ...)     standardGeneric("plot.cghSeg"))
setGeneric("segmented",     function(object)        standardGeneric("segmented"))
setGeneric("segmented<-",   function(object, value) standardGeneric("segmented<-"))

### cghCall

setGeneric("plot.cghCall",  function(x, y, ...)     standardGeneric("plot.cghCall"))
#setGeneric("summaryPlot",   function(x, y, ...)     standardGeneric("summaryPlot"))
setGeneric("calls",         function(object)        standardGeneric("calls"))
setGeneric("calls<-",       function(object, value) standardGeneric("calls<-"))
setGeneric("probdloss",      function(object) standardGeneric("probdloss"))
setGeneric("probdloss<-",    function(object, value) standardGeneric("probdloss<-"))
setGeneric("probloss",      function(object)        standardGeneric("probloss"))
setGeneric("probloss<-",    function(object, value) standardGeneric("probloss<-"))
setGeneric("probnorm",      function(object)        standardGeneric("probnorm"))
setGeneric("probnorm<-",    function(object, value) standardGeneric("probnorm<-"))
setGeneric("probgain",      function(object)        standardGeneric("probgain"))
setGeneric("probgain<-",    function(object, value) standardGeneric("probgain<-"))
setGeneric("probamp",       function(object)        standardGeneric("probamp"))
setGeneric("probamp<-",     function(object, value) standardGeneric("probamp<-"))

### cghRegions

setGeneric("plot.cghRegions",   function(x, y, ...)     standardGeneric("plot.cghRegions"))
setGeneric("frequencyPlot",     function(x, y, ...)     standardGeneric("frequencyPlot"))
setGeneric("nclone",            function(object)        standardGeneric("nclone"))
setGeneric("avedist",           function(object)        standardGeneric("avedist"))
setGeneric("regions",           function(object)        standardGeneric("regions"))
setGeneric("regions<-",         function(object, value) standardGeneric("regions<-"))

#setGeneric("plot",          function(x, y, ...) standardGeneric("plot"), useAsDefault=plot)
setGeneric("plot",          function(x, y, ...) standardGeneric("plot"), useAsDefault=plot)
