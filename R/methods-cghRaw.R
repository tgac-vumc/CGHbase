setMethod("initialize", "cghRaw",
        function(.Object,
                assayData       = assayDataNew(copynumber=copynumber, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHbase:::.makeEmptyFeatureData(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                copynumber      = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setValidity("cghRaw", function(object) {
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "cghRaw"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("copynumber")))
    msg <- Biobase:::validMsg(msg, CGHbase:::.featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End")))
    if (is.null(msg)) TRUE else msg
})

setMethod("plot", signature(x="cghRaw", y="missing"),
function (x, y, dotres=10, ylimit=c(-2,5), ylab=expression(log[2]~ratio), build="GRCh37",... )
{
    for (i in 1:ncol(x)) {
        cat("Plotting sample", sampleNames(x)[i], "\n")
        chrom           <- chromosomes(x)
        data            <- data.frame(chrom, bpstart(x), copynumber(x)[,i])
        colnames(data)  <- c("chromosome", "position", "ratio")
        pos             <- bpstart(x)
        uni.chrom <- unique(chrom)
        chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
        chrom.ends <- integer()
        cumul <- 0
        for (j in uni.chrom) {
            pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
            cumul <- cumul + chrom.lengths[as.character(j)]
            chrom.ends <- c(chrom.ends, cumul)
        }
        names(chrom.ends) <- names(chrom.lengths)
        nclone <- length(chrom)
        whichtoplot <- seq(1,nclone,by=dotres) #added 15/06/2009
        plot(pos[whichtoplot], data[whichtoplot,3], cex=.1, main=sampleNames(x)[i], ylab=ylab, xlab="chromosomes", ylim=ylimit, xaxt="n", xaxs="i")
        if (dotres != 1)
            mtext(paste('Plot resolution: ', 100/dotres, '%', sep=''), side=3, line=0)
        abline(h=0)
        if (length(chrom.ends) > 1)
            for (j in names(chrom.ends)[-length(chrom.ends)])
                abline(v=chrom.ends[j], lty='dashed')
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1, at=ax, labels=uni.chrom, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1)
        amps <- data[,3]
        amps[amps>=5] <- 5.15
        amps[amps<5] <- NA
        points(pos, amps, pch=24, col='blue', bg='blue', cex=0.5)
        dels <- data[,3]
        dels[dels<=-2] <- -2.15
        dels[dels>-2] <- NA
        points(pos, dels, pch=25, col='red', bg='red', cex=0.5)
        ### MAD
        windowsize <- 50
        x1 <- copynumber(x)[chromosomes(x) < 23,i]
        elc <- length(x1)-200
        
        if(elc>=1000){
        seqs <- seq(1,elc, by=floor(elc/100))
        mad.value <- round(median(sapply(seqs,function(wh) mad(x1[wh:(wh+windowsize)], na.rm=TRUE))), digits=2)
        mtext(paste('MAD =', mad.value), side=3, line=0, adj=1)
        }
        ### number of data points
        str <- paste(round(nclone / 1000), 'k x ', sep='')
        probe <- median(bpend(x)-bpstart(x)+1)
        if (probe < 1000) {
            str <- paste(str, probe, ' bp', sep='')
        } else {
            str <- paste(str, round(probe / 1000), ' kbp', sep='')
        }
        mtext(str, side=3, line=0, adj=0)
    }
})



setMethod("chromosomes", "cghRaw", function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart",     "cghRaw", function(object) pData(featureData(object))[,"Start"])
setMethod("bpend",       "cghRaw", function(object) pData(featureData(object))[,"End"])

setMethod("copynumber", signature(object="cghRaw"),
        function(object) Biobase:::assayDataElement(object, "copynumber"))
        
setReplaceMethod("copynumber", signature(object="cghRaw", value="matrix"),
                function(object, value) assayDataElementReplace(object, "copynumber", value))    
