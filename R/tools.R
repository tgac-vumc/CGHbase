cghRaw <- function(input) {
    if (class(input) == "character") {
        input   <- read.table(input, header=T, sep="\t", fill=T, quote="")
    } else if (class(input) == "data.frame") {
        input   <- input
    } else {
        cat("Input should be either a data.frame or a character string containing a filename.\nPlease read the documentation file for the required format.\n")
        input   <- NULL
    }
    if (!is.null(input)) {
        copynumber  <- as.matrix(input[,5:ncol(input)])
        rownames(copynumber) <- input[,1]
        annotation  <- data.frame(Chromosome=input[,2], Start=input[,3], End=input[,4], row.names=input[,1])
        metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))    
        dimLabels   <- c("featureNames", "featureColumns")
        annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)   
        result      <- new("cghRaw", copynumber=copynumber, featureData=annotation)
    }
}

.setFeatureData <- function(object) {
    ### To be developed when a dataset with annotation is available
    lib <- annotation(object)
    require(package=lib, character.only=TRUE)
    
}
