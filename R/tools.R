make_cghRaw <- function(input) {
    if (class(input) == "character") input  <- read.table(input, header=T, sep="\t", fill=T, quote="")
    if (class(input[,2]) == 'factor')
        input[,2] <- as.character(input[,2])
    if (class(input[,2]) == 'character') {
        input[,2] <- sub('^chr', '', input[,2])
        input[input[,2] == 'X', 2] <- '23'
        input[input[,2] == 'Y', 2] <- '24'
        input[input[,2] == 'MT', 2] <- '25'
        input[,2] <- as.integer(input[,2])
    }
    input <- input[order(input[,2], input[,3]),]
    copynumber  <- as.matrix(input[,5:ncol(input)])
    rownames(copynumber) <- input[,1]
    if (ncol(copynumber) == 1)
      colnames(copynumber) <- colnames(input)[5]
    annotation  <- data.frame(Chromosome=input[,2], Start=input[,3], End=input[,4], row.names=input[,1])
    metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))    
    dimLabels   <- c("featureNames", "featureColumns")
    annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)   
    result      <- new("cghRaw", copynumber=copynumber, featureData=annotation)
}

.setFeatureData <- function(object) {
    ### To be developed when a dataset with annotation is available
    lib <- annotation(object)
    require(package=lib, character.only=TRUE)
    
}
