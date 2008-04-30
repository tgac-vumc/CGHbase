.featureDataRequiredColumns <- function(featureData, columns) {
    msg     <- NULL
    absent  <- columns[!(columns %in% rownames(varMetadata(featureData)))]
    if (length(absent) != 0) {
        msg <- paste(msg, paste("missing columns' ", absent ,"' in featureData" , sep = "", collapse = "\n\t"), sep="\n")
    }
    if (is.null(msg)) TRUE else msg
}

.makeEmptyFeatureData <- function(object) {
    dims        <- Biobase:::assayDataDims(object)
    n           <- dims[1,1]
    features    <-         
    if (is(object, "environment")) ls(object)
    else names(object)
    nms         <- rownames(object[[features[[1]]]])
    data        <- data.frame(Chromosome=numeric(n), Start=numeric(n), End=numeric(n), row.names=nms)
    dimLabels   <- c("featureNames", "featureColumns")
    metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))
    new("AnnotatedDataFrame", data=data, dimLabels=dimLabels, varMetadata=metadata)                  
}

.makeEmptyFeatureDataForRegions <- function(object) {
    dims        <- Biobase:::assayDataDims(object)
    n           <- dims[1,1]
    features    <-         
    if (is(object, "environment")) ls(object)
    else names(object)
    nms         <- rownames(object[[features[[1]]]])
    data        <- data.frame(Chromosome=numeric(n), Start=numeric(n), End=numeric(n), Nclone=numeric(n), AveDist=numeric(n), row.names=nms)
    dimLabels   <- c("featureNames", "featureColumns")
    metadata    <- data.frame(  labelDescription=c("Chromosomal position", 
                                                    "Basepair position start", 
                                                    "Basepair position end", 
                                                    "Number of clones in region", 
                                                    "Average distance"), 
                                row.names=c("Chromosome", 
                                            "Start", 
                                            "End",
                                            "Nclone",
                                            "AveDist")
                                )
    new("AnnotatedDataFrame", data=data, dimLabels=dimLabels, varMetadata=metadata)                  
}

.makeSegments <- function(data) {
    if (class(data) != "matrix") {
        cat("Wrong data class for input to function CGHbase:::.makeSegments in file private.R\n")
    }
    previous    <- 2000
    values      <- c()
    start       <- c()
    end         <- c()
    for (i in 1:nrow(data)) {
        if (!all(data[i,] == previous)) {
            start   <- c(start, i)
            last    <- i - 1
            if (last > 0) end <- c(end, last)
            values  <- c(values, data[i,2])
        }
        previous    <- data[i,]
    }
    end     <- c(end, nrow(data))
    result  <- cbind(values, start, end)
    result
}
