setMethod("initialize", "cghRegions",
        function(.Object,
                assayData       = assayDataNew(regions=regions, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHbase:::.makeEmptyFeatureDataForRegions(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                regions         = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setMethod("plot.cghRegions", signature(x="cghRegions", y="missing"),
function (x, y, ... )
{
    regions <- x
    
    prlossgain <- function(aber, sign, bound, nsam) {
        numga   <- (length(sign[sign == aber]))/nsam
        bound1  <- c(0, bound, 1)
        for (i in 1:(length(bound)+1)) {
            if (numga >= bound1[i] & numga <= bound1[i+1]) clas <- i
        }
        return(clas)
    }    
    
    nsam        <- ncol(regions)
    boundary    <- c(0.1,0.3,0.5)
    prlosscl    <- apply(regions(regions), 1, prlossgain, aber=-1, bound=boundary, nsam=nsam)
    prgaincl    <- apply(regions(regions), 1, prlossgain, aber=1, bound=boundary, nsam=nsam)
    palloss     <- marray::maPalette(low = "black", high = "red", k=4)
    palgain     <- marray::maPalette(low = "black", high = "green", k=4)
    colloss     <- sapply(prlosscl, function(a){return(palloss[a])})
    colgain     <- sapply(prgaincl, function(a){return(palgain[a])})    
    
    xall    <- c()
    x1all   <- c()
    x2all   <- c()
    yall    <- c()
    y1all   <- c()
    for(i in unique(chromosomes(regions))) {
        datachr <- regions[chromosomes(regions) == i,]
        nr      <- nrow(datachr)
        x1      <- bpstart(datachr)
        x2      <- bpend(datachr)
        y1      <- c()
        for(j in 1:nr) {
            y1 <- c(y1,23-i+0.16*(-1)^(j+1))
        }
        x       <- c(x1,x2)
        y       <- c(y1,y1)
        xall    <- c(xall,x)
        x1all   <- c(x1all,x1)
        x2all   <- c(x2all,x2)
        yall    <- c(yall,y)
        y1all   <- c(y1all,y1)
    } 

    layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE),heights=c(8,1))
    par(mar=c(2.5,3,1,1.5))
    plot(xall, yall, cex=0.1, yaxt="n", ylab="", xlab="")
    delta       <- 0.17
    yallmin     <- yall - delta
    y1allmin    <- y1all - delta
    points(xall, yallmin, cex=0.1)
    segments(x1all, y1all, x2all, y1all, col=colgain, lwd=5, yaxt="n")  
    segments(x1all, y1allmin, x2all, y1allmin, col=colloss, lwd=5, yaxt="n") 
    axis(2, at=1:22, labels=22:1, las=1, cex.axis=1, cex.lab=1)
    marray::maColorBar(1:4, col=palloss, horizontal=TRUE, k=4,yaxt="n")
    marray::maColorBar(1:4, col=palgain, horizontal=TRUE, k=4,yaxt="n")
})

setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
function (x, y, ... )
{  
    regions     <- x
    
    freqlossgain <- function(aber, sign, nsam) {
        numga <- (length(sign[sign == aber]))/nsam
        return(numga)
    }
    
    nsam        <- ncol(regions)
    proplosscl  <- apply(regions(regions), 1, freqlossgain, aber=-1, nsam=nsam)
    propgaincl  <- apply(regions(regions), 1, freqlossgain, aber=1, nsam=nsam)
    
    par(mar=c(5, 4, 4, 4) + 0.2)
    
    propnormcl  <- 1-(proplosscl+propgaincl)
    probsdraw   <- cbind(proplosscl, propnormcl, propgaincl)
    chr         <- chromosomes(regions)
    bp          <- cbind(bpstart(regions), bpend(regions))
    bplength    <- bp[,2]-bp[,1]
    totalwidth  <- cumsum(as.numeric(bplength))
    tlength     <- totalwidth[length(totalwidth)]
    bplsc       <- bplength/tlength
    twsc        <- totalwidth/tlength
    chrbor      <- cbind(chr,twsc)
    x2all       <- c()
    
    for(i in unique(chr)) {
        datachr <- chrbor[chrbor[,1]==i,]
        nr      <- nrow(datachr)
        x2all   <- c(x2all,if(!is.null(nr)){datachr[nr,2]} else {datachr[2]})
    } 

    barplot(t(probsdraw), width=bplsc, space=0, border=F, col=c("gray", "white", "black"), las=1, cex.axis=1, cex.lab=1, xaxt="n", xaxs="i")
    for (iii in 1:length(unique(chr))) {segments(x2all[iii], -5, x2all[iii], 5, lty=2, col=gray(0.5))}
    lim   <- par("usr")
    lim[3:4] <- c(-5,5)
    par(usr=lim)
    box()
    ax <- (x2all + c(0, x2all)[-(length(x2all)+1)])/2
    axis(side=1, at=ax, labels=unique(chr), cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1) # bottom axis  
})

setValidity("cghRegions", function(object) {
    msg <- NULL
    msg <- Biobase:::validMsg(msg, Biobase:::isValidVersion(object, "cghRegions"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("regions")))
    msg <- Biobase:::validMsg(msg, CGHbase:::.featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End", "Nclone", "AveDist")))
    if (is.null(msg)) TRUE else msg
})

setMethod("chromosomes",    "cghRegions", function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart",        "cghRegions", function(object) pData(featureData(object))[,"Start"])
setMethod("bpend",          "cghRegions", function(object) pData(featureData(object))[,"End"])
setMethod("nclone",         "cghRegions", function(object) pData(featureData(object))[,"Nclone"])
setMethod("avedist",        "cghRegions", function(object) pData(featureData(object))[,"AveDist"])

setMethod("regions", signature(object="cghRegions"),
        function(object) Biobase:::assayDataElement(object, "regions"))
        
setReplaceMethod("regions", signature(object="cghRegions", value="matrix"),
                function(object, value) assayDataElementReplace(object, "regions", value))  
