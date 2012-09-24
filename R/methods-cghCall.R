setMethod("initialize", "cghCall",
        function(.Object,
                assayData       = assayDataNew(copynumber=copynumber, segmented=segmented, calls=calls, probloss=probloss, probnorm=probnorm, probgain=probgain, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHbase:::.makeEmptyFeatureData(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                copynumber      = new("matrix"),
                segmented       = new("matrix"),
                calls           = new("matrix"),
                probloss        = new("matrix"),
                probnorm        = new("matrix"),
                probgain        = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setValidity("cghCall", function(object) {
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "cghCall"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("copynumber", 
                                                                                        "segmented", 
                                                                                        "calls", 
                                                                                        "probloss", 
                                                                                        "probnorm", 
                                                                                        "probgain"
                                                                                        )))
    msg <- Biobase:::validMsg(msg, .featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End")))
    if (is.null(msg)) TRUE else msg
})

#probdloss <- function(object) Biobase:::assayDataElement(object, "probdloss")

setMethod("plot", signature(x="cghCall", y="missing"),
function (x, y, dotres=10, ylimit=c(-5,5), ylab=expression(log[2]~ratio), gaincol='green', losscol='red',ampcol="darkgreen",dlcol="darkred", misscol=NA, build='GRCh37',... )
{
#x<-calls2; dotres=100; ylimit=c(-5,5); ylab=expression(log[2]~ratio); gaincol='green'; losscol='red'; ampcol="darkgreen";dlcol="darkred"; misscol=NA; build='GRCh37'
    calls           <- calls(x)
    nsamples        <- ncol(x)
    nclass <-3
    if (!is.null(probamp(x))) nclass <- nclass+1 
    if (!is.null(probdloss(x))) nclass <- nclass+1 
    chrom           <- chromosomes(x)
    pos             <- bpstart(x)
    pos2            <- bpend(x)
    uni.chrom <- unique(chrom)
    chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
    chrom.ends <- integer()
    cumul <- 0
    for (j in uni.chrom) {
        pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
        pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
        cumul <- cumul + chrom.lengths[as.character(j)]
        chrom.ends <- c(chrom.ends, cumul)
    }
    names(chrom.ends) <- names(chrom.lengths)
    nclone          <- length(chrom)
    chrom.labels    <- unique(chrom)
    
   

     for (i in 1:ncol(x)) {
     #i<-1
        cat("Plotting sample", sampleNames(x)[i], "\n")
        inc <- 0
        genomdat        <- copynumber(x)[,i]
        if(nclass==3) probsdraw       <- cbind(probloss(x)[,i], probnorm(x)[,i], probgain(x)[,i])
        if(nclass==4) probsdraw       <- cbind(probloss(x)[,i], probnorm(x)[,i], probgain(x)[,i]+probamp(x)[,i])
        if(nclass==5) {
            probsdraw       <- cbind(probloss(x)[,i]+probdloss(x)[,i], probnorm(x)[,i], probgain(x)[,i]+probamp(x)[,i])
            inc <- 1
            }
        lt              <- 0
        ltdl            <- 0
        
        if (nclass>=4) {
            ticksamp    <- which(probamp(x)[,i] >= 0.5)
            lt          <- length(ticksamp)        
        }
        
        if (nclass==5) {
            ticksdl    <- which(probdloss(x)[,i] >= 0.5)
            ltdl          <- length(ticksdl)
            #print(ltdl)        
        }
        
        segment         <- .makeSegments(segmented(x)[,i],chromosomes(x))
        
        widths          <- segment[,3] - segment[,2] + 1
        par(mar=c(5, 4, 4, 4) + 0.2)
        
        ### Plot the probability bars
        plot(NA, xlim=c(0, max(pos2)), ylim=c(0,1), xlab=NA, ylab=NA, las=1, xaxs='i', xaxt='n', yaxs='i')
        if (!is.na(misscol)) {
            rect(0, -1, max(pos2), 1, col=misscol, border=NA)
            rect(pos, -1, pos2, 1, col='white', border=NA)
        }
        rect(pos[segment[,2]], 0, pos2[segment[,3]], probsdraw[segment[,2],1], col=losscol, border=losscol)
        rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-probsdraw[segment[,2],3], col=gaincol, border=gaincol)
        
        lim <- par("usr")
        lim[3:4] <- ylimit
        par(usr=lim)
        dticks <- seq(ylimit[1], ylimit[2], by=1)
        axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1, cex.lab=1)
 
        if (lt > 0) {
            axis(3,at=pos[ticksamp], labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1, cex.axis=1, cex.lab=1)
        }  
        
        if (ltdl > 0) {
            axis(3,at=pos[ticksdl], labels=FALSE, col=dlcol, col.axis="black", srt=270, las=1, cex.axis=1, cex.lab=1)
        }   
        
               
        box()   
        abline(h=0)
        
        ### Add axis labels
        mtext(ylab, side=4, line=3, srt=270)
        mtext("probability", side=2, line=3, srt=270)
        mtext('chromosomes', side=1, line=3)
        
        #### add vert lines at chromosome ends
        if (length(chrom.ends) > 1)
            for (j in names(chrom.ends)[-length(chrom.ends)])
                abline(v=chrom.ends[j], lty='dashed')

        title(sampleNames(x)[i])
        if (dotres != 1)
            mtext(paste('Plot resolution: ', 100/dotres, '%', sep=''), side=3, line=0)
        
        ### Add log2ratios
        if (dotres>0) {
            whichtoplot <- seq(1,nclone,by=dotres) #added 15/06/2009
            points(pos[whichtoplot],genomdat[whichtoplot],cex=.1)
        }
        
        ### X-axis with chromosome labels
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
            
        ### segment means
        for (jjj in (1:nrow(segment)))
            segments(pos[segment[jjj,2]], segment[jjj,1], pos[segment[jjj,3]], segment[jjj,1], col="chocolate", lwd=3)        

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

frequencyPlotCalls <- 
function(x, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  calls <- calls(x)
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
}



summaryPlot <- 
function (x, main='Summary Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  
  nclass <-3
  if (!is.null(probamp(x))) nclass <- nclass+1 
  if (!is.null(probdloss(x))) nclass <- nclass+1 
  
  chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  
  if(nclass==3) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))}
  if(nclass==4) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}
  if(nclass==5) {loss.freq <- rowMeans(probloss(x))+rowMeans(probdloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}
  
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='mean probability', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
}


setMethod("chromosomes", "cghCall", function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart",     "cghCall", function(object) pData(featureData(object))[,"Start"])
setMethod("bpend",       "cghCall", function(object) pData(featureData(object))[,"End"])

setMethod("copynumber", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "copynumber"))
        
setReplaceMethod("copynumber", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "copynumber", value))   

setMethod("segmented", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "segmented"))
        
setReplaceMethod("segmented", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "segmented", value))  
                
setMethod("calls", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "calls"))
        
setReplaceMethod("calls", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "calls", value))                  

setMethod("probdloss", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probdloss"))
        
setReplaceMethod("probdloss", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probdloss", value))

setMethod("probloss", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probloss"))
        
setReplaceMethod("probloss", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probloss", value))
                
setMethod("probnorm", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probnorm"))
        
setReplaceMethod("probnorm", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probnorm", value))
                
setMethod("probgain", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probgain"))
        
setReplaceMethod("probgain", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probgain", value))    
                
setMethod("probamp", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probamp"))
        
setReplaceMethod("probamp", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probamp", value))  
