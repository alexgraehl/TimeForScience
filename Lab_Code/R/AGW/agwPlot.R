
## =============================================================================
## This file contains various miscellaneous plotting functions. They should be general-purpose.
## =============================================================================

if (!exists("print.agw")) {
     if (file.exists("/work/Common/Code/R_Binf_Core/Utility/agwUtil.R")) {
          source("/work/Common/Code/R_Binf_Core/Utility/agwUtil.R")
     } else if (file.exists("~/TimeForScience/Lab_Code/R/AGW/agwUtil.R")) {
          source("~/TimeForScience/Lab_Code/R/AGW/agwUtil.R")
     } else if (file.exists("./agwUtil.R")) {
          source("./agwUtil.R")
     } else {
          stop("Error -- can't find Alex's R utilities on the filesystem.")
     }                      
}


## Plotting-related functions.

### ===============================================================================
## Just a plain rectangle with a bunch of (wrapped) text in it.
## wraplen: the number of characters before we wrap text
## leftMartin: the left margin. Ranges from 0-1. 0 = no margin, 1 = figure is only margin, and no text can be printed
## topMargin: the top start location. Ranges from 0-1. 0 means "start at the top, no margin" and 1.0 means "don't show the text, it's all margin"
### ===============================================================================
textDescriptionPlotAGW <- function(text, wraplen=60, leftMargin=0.15, topMargin=0.01, cex=1) {
     prevMar <- par()$mar  ;  par(mar=c(0,0,0,0)) ## no margins for this plot...
     wordwrap <- function(str, len) {
          modStr <- (strsplit(str, "\n"))[[1]]
          return(paste(strwrap(modStr, width=len), collapse="\n"))
     }
     plot(c, axes=F, xlab=NA, ylab=NA, type='n', frame.plot=F)
     text(x=leftMargin, y=(1-topMargin), adj=c(0,1), cex=cex, labels=wordwrap(text, len=wraplen))
     par(mar=prevMar) ## restore the old margins...
}





## =================================================================
## Colors -- gradient, like heat.colors
## =================================================================
colors.agw <- function(n = 12, type="blueblackyellow", reverse=FALSE) {
     # Returns a color gradient, much like "heat.colors(...)" Several options for "type" are available.
     stopifnot(is.numeric(n))
     stopifnot(n >= 1)
     
     totalColors <- n
     colRange01 <- (0:(totalColors-1))/(totalColors-1)
     range1 <- (0:(ceiling(totalColors / 2)-1))/(ceiling(totalColors/2)-1) ## For two-color gradients (the left half of the colors)
     range2 <- (1:(floor(totalColors / 2)))/(floor(totalColors/2))         ## For two-color gradients (the right half of the colors)

     type <- tolower(type)

     if (type == "greenwhitered") { ## green -> white -> red
          col1 <- rev(hsv(h=rev(0.3+0.20*range1), s=range1, v=(1.0 - 0.2*range1)))
          col2 <- (hsv(h=rev(0.0+0.15*range2), s=range2, v=(1.0 - 0.15*range2)))
          col <- c(col1, col2)
     } else if (type == "greenblackred") { ## green -> black -> red
          col1 <- rev(hsv(h=rev(0.3+0.20*range1), s=rev(1.0-0.2*range1), v=range1))
          col2 <- (hsv(h=rev(0.0+0.10*range2), s=rev(1.0-0.2*range2), v=range2))
          col <- c(col1, col2)
     } else if (type == "blueblackyellow") { ## blue -> black -> yellow
          col1 <- rev(hsv(h=rev(0.6+0.20*range1), s=rev(1.0-0.2*range1), v=range1))
          col2 <- (hsv(h=rev(0.15+0.10*range2), s=rev(1.0-0.2*range2), v=range2))
          col <- c(col1, col2)
     } else if (type == "blueblackyellow2") { ## blue -> black -> yellow, hand-picked
          stopifnot(n == 11)
          col <- c("#00C8FF", "#00AAF5", "#0082D7", "#2464A8", "#004064", "black", "#646400","#919114","#B6B61E","#D7D728","#FFFF00")
     } else if (grepl("^gr[ea]y", type)) { col <- gray(colRange01) }
     else if (grepl("^brown", type))   { col <- rev(hsv(h=rev(0.2*colRange01), s=colRange01, v=(1.0 - 0.7*colRange01))) }
     else if (grepl("^sepia", type))   { col <- rev(hsv(h=rev(0.3*colRange01), s=colRange01, v=rev(colRange01))) }
     else if (grepl("^heat", type))    { col <- heat.colors(n) }
     else {
          print(paste("An unrecognized color gradient type (<", type, ">) was passed into colorGradientAGW:", sep=''))
          print(type)
          stopifnot(paste("Color type is not recognized. Try something like \"gray.colors\"") == 999)
     }
     if (reverse) { col <- rev(col) }
     return(col)
}




## =================================================================
## pdf.for.heatmap.agw:
## * This is used INSTEAD of a call to pdf
## * Automatically calculates the proper PDF dimentions for a heatmap.agw plot.
## * This way the heatmaps can be scaled to have cells that are semi-approximately the right size.
##   Otherwise, heatmaps with lots of rows tend to have really cramped rows, and heatmaps with few rows
##   have gigantic super-tall cells.
## * You can pass in either the matrix that you are about to plot, OR the numRows you will plot.
## =================================================================
pdf.for.heatmap.agw <- function(file=file, mat=NULL, numRows=NULL, width="should not be specified by the user!!", height="should not be specified by the user!!", ...) {
     if (!is.null( mat )) {
          amount <- nrow( mat )
          assert.agw(is.null(numRows), paste("You cannot specify both a matrix AND a number of rows!! Input one OR the other! numRows was specified to be: ", numRows, ".", sep=''))
     } else if (!is.null(numRows)) {
          amount <- numRows
     }
     assert.agw(missing(width), "Hey! You should not specify width in pdf.for.heatmap.agw! It is always set to the same value.")
     assert.agw(missing(height), "Hey! You cannot manually specify height in pdf.for.heatmap.agw. It is auto-calculated!")
     assert.agw(is.numeric(amount), "Uh oh! Pdf.for.heatmap.agw is broken!")

     ALEX_HEATMAP_WIDTH_INCHES <- 15
     MAX_HEIGHT <- 80
     MIN_HEIGHT <- 12
     pdf(file=file,
         , width=ALEX_HEATMAP_WIDTH_INCHES
         , height=min(MAX_HEIGHT, max(6+amount*0.25, 12))
         , ...)
}

## =================================================================
## Alex's Heatmap-and-histogram plot
## For a PDF, this heatmap.agw must be at least 12 inches tall for the heatmap AND the histogram to both fit.
## It can be 8 or more inches wide and look OK.
## =================================================================
heatmap.agw <- function(mmm, breaks=12, labRow=NULL, labCol=NULL, col=NULL, colorStyle=NULL, main, title="", cexRow=NULL, cexCol=1.5, maxNumLabels=1000, col.names=NULL, row.names=NULL, cluster.rows=FALSE) {
     ## M: a matrix to plot
     ## Breaks: the number of histogram breaks, used for the color scheme
     ## maxNumLabels: do not print labels if there are more than this many labels ***with actual non-blank content***

     print.agw("heatmap.agw: Now generating a \"heatmap.agw\" figure. If you get a \"figure region too large\" error,")
     print.agw("             that means your PDF/PNG wasn't big enough to fit the heatmap. This can be solved by using")
     print.agw("             a bigger pdf width or by using \"pdf.for.heatmap.agw\" to auto-compute the bounds.")
     
     if (is.logical(cluster.rows) && cluster.rows) {
          ## If cluster.rows is true, then we will CLUSTER the rows, kind of like how regular built-in "heatmap"
          ## does it. Check the source to "heatmap" to get a sort of general idea how this works.
          hcc <- stats::hclust(stats::dist(mmm, method="euclidean")
                               , method="complete")
          dcc <- as.dendrogram(hcc)
          reordering <- order.dendrogram(dcc)
          mmm <- mmm[reordering, , drop=F] ## REORDER THE INPUT MATRIX BASED ON THE CLUSTERING
               assert.agw(is.null(labRow) && is.null(row.names), "Uh oh! You cannot specify that you want the matrix to be re-clustered AND ALSO specify row names. This is because once you recluster, the row names will not be what you probably expect! i.e., the row names move around!")
     }
     
     ## col: the SPECIFIC list of colors to pass in. Must be equal in length to (breaks - 1)
     ## colorStyle: OR you can specify the colors as a style. This is one of the strings accepted by "colors.agw"--for example, "sepia" or "gray" or "blueblackyellow"

     if (is.vector(mmm)) {
          mmm <- as.matrix(mmm)
     }
     stopifnot(is.matrix(mmm))
     stopifnot(nrow(mmm) >= 1)
     stopifnot(ncol(mmm) >= 1)
     
     # Generates a three-row multi-part figure.
     # Top part: the caption (title)
     # Middle part: the histogram and distribution key (probably should be made optional, but it's required for now)
     # Bottom part: the heatmap
     
     if (!missing(col.names) && !is.null(col.names)) { labCol <- col.names } ## "col.names" is just an alias for "labCol"
     if (!missing(row.names) && !is.null(col.names)) { labRow <- row.names } ## "row.names" is just an alias for "labRow"
     
     if (is.null(labRow)) { labRow = rownames(mmm) }
     if (is.null(labCol)) { labCol = colnames(mmm) }


     stopifnot(breaks >= 3)

     scale01 <- function(x, low = min(x), high = max(x)) { return((x - low)/(high - low)) }

     layout(matrix(c(1,2,3), byrow=T, ncol=1), heights=c(lcm(6),lcm(12),1) ) # <-- we plot THREE things, stacked vertically

     # ==========================================
     # ======================================

     min.raw <- min(mmm, na.rm=TRUE)
     max.raw <- max(mmm, na.rm=TRUE)

     mean.raw <- mean(mmm, na.rm=TRUE)
     median.raw <- median(mmm, na.rm=TRUE)

     ## This is the FIRST of three "layout" sub-items (which are stacked vertically). It's a text description of what is being plotted.
     textDescriptionPlotAGW(paste("Heatmap with ", length(mmm)
                                  , " values, in ", nrow(mmm)
                                  , " rows and "
                                  , ncol(mmm), " columns."
                                  , " Mean = ", format(mean.raw, digits=3, nsmall=1)
                                  , ", median = ", format(median.raw, digits=3, nsmall=1)
                                  , ". ", title
                                  , sep=''), wraplen=120, leftMargin=0.05, topMargin=0.05, cex=1.4)

     # ======================================
     # ==========================================
     # ==========================================
     # ======================================

     ## This is the SECOND of three layout things. It's a histogram
     
     par(mar=c("bottom"=3, "left"=9.5, "top"=5, "right"=15.5))

     if (missing(main) || is.null(main)) {
          main <- "Heatmap"
     }
     main <- paste(main, "\n(", nrow(mmm), " rows by ", ncol(mmm), " columns)", sep='')

     if (missing(breaks) || is.null(breaks) || (length(breaks) == 1)) {
          ## If breaks was not specified, OR it was a length-one scalar
          breaks <- seq(min.raw, max.raw, length.out=breaks)
     }

     ## Alternative color scale chosen by Alex, blue to yellow (black in the middle):
     ##, col=c("#00C8FF", "#00AAF5", "#0082D7", "#2464A8", "#004064", "black", "#646400","#919114","#B6B61E","#D7D728","#FFFF00") ## blue to yellow

     if (!missing(colorStyle) && !is.null(colorStyle)) {
          assert.agw(length(colorStyle) == 1, "colorStyle needs to be a string like blueblackyellow. The user can pass in a string here to automagically pick the color scheme. Or they can specify it manually with \"col\".")
          ## ColorStyle is a CHARACTER vector
          col <- colors.agw(n = length(breaks)-1, type=colorStyle)
     }

     if (missing(cexCol) || is.null(cexCol)) {
          cexCol = 1.5
          if (nrow(mmm) > 12) {
               cexCol=1.0
          }
          if (nrow(mmm) > 25) {
               cexCol=0.8
          }
     }
     
     if (missing(col) || is.null(col)) {
          ## If the color was not specified...
          totalColors <- length(breaks)-1
          if (totalColors == 2) {
               col <- c("#000066", "#FFFF99")
          } else if (totalColors >= 10) {
               col <- c("#330033", "#440044", "#550055", "#770044", "darkred", heat.colors(totalColors-6))
          } else {
               col <- c("#330033", "#770044", heat.colors(totalColors-2))
          }
          col <- col[1:totalColors]
     }

     ## Draw the KEY / LEGEND histogram
     
     ## Draw the background for the "legend" histogram
     legendBack <- seq(min.raw, max.raw, length=length(col)) ## Heatmap Histogram!
     image(z = matrix(legendBack, ncol = 1), col=col, breaks=breaks, xaxt="n", yaxt="n") ## This is the histogram / distribution background that goes at the top of the heatmap
     
     par(usr = c(0, 1, 0, 1))
     lv <- pretty(breaks)
     xv <- scale01(as.numeric(lv), min.raw, max.raw)
     axis(1, at=xv, labels=lv)

     h <- hist(mmm, plot = FALSE, breaks=breaks)
     hx <- scale01(breaks, min.raw, max.raw)
     hy <- c(h$counts, h$counts[length(h$counts)])

     #print(mean.raw)
     #abline(v=scale01(mean.raw, min.raw, max.raw), lwd=4, col="black")

     MEDIAN_LINE_WIDTH    <- 4
     HISTOGRAM_LINE_WIDTH <- 6
     #print(median.raw)
     abline(v=scale01(median.raw, min.raw, max.raw), lwd=2*MEDIAN_LINE_WIDTH, lty="solid", col="gray")
     abline(v=scale01(median.raw, min.raw, max.raw), lwd=MEDIAN_LINE_WIDTH, lty="dashed", col="black")

     Y_SCALE <- 0.95
     lines(hx, hy/max(hy) * Y_SCALE, lwd=2*HISTOGRAM_LINE_WIDTH, type = "s", col = "white") ## Draw the actual histogram as a squiggly line
     lines(hx, hy/max(hy) * Y_SCALE, lwd=HISTOGRAM_LINE_WIDTH, type = "s", col = "black")

     axis(2, at=pretty(hy)/max(hy) * Y_SCALE, pretty(hy))
     par(cex.main=1.7)
     title(paste("Color Key and Histogram of the Distribution of Heatmap Values\nDashed line = median value (", format(median.raw, digits=3, nsmall=1), ")", sep=''))
     mtext(side=2, paste("Count (out of ", length(mmm) ,")", sep=''), line=3)
     box(lwd=1)
     ## Done drawing the "legend" histogram of values in the heatmap
     
     # ======================================
     # ==========================================
     # ==========================================
     # ======================================

     ## This is the THIRD thing being plotted, and is the real heatmap. If you want a heatmap that doesn't use layout/mfrow/etc., you can just copy out this code here.
     ## This is the "meat" of the heatmap generation.
     
     par(mar=c("bottom"=25, "left"=8, "top"=8, "right"=20), cex.main=2.2)
     par(cex = 0.5) ## Expansion factor
     #image(m, axes=F, main=main, col=col)
     image(t(mmm)[, nrow(mmm):1, drop=F], axes=F, main=main, col=col) ## <-- transpose AND flip, to rotate the correct way
     # Notice that ‘image’ interprets the matrix as a table of
     # ‘f(x[i], y[j])’ values, so that the x axis corresponds to row
     # number and the y axis to column number, with column 1 at the
     # bottom, i.e. a 90 degree counter-clockwise rotation of the
     # conventional printed layout of a matrix.

     numNonBlankRows <- 0
     if (!is.null(labRow)) {
          numNonBlankRows <- sum(!is.na(labRow) & (nchar(labRow) > 0)) ## Count the number of NON-BLANK rows only!
     }
     
     if (is.null(cexRow)) {
          rowNumberToUseForSizeCalculation <- nrow(mmm) ## Make the labels tiny enough to individually specify a single row
          cexRow = 1.0
          rowsAtWhichCexIsMin <- 800 ## The number of rows at which the text is the tiniest (or more rows than this)
          rowsAtWhichCexIsMax <- 100 ## The number of rows at which the text is the biggest (or fewer rows than this)
          stopifnot(rowsAtWhichCexIsMin > rowsAtWhichCexIsMax)
          maxRowCex <- 1.00
          minRowCex <- 0.10

          if (rowNumberToUseForSizeCalculation >= rowsAtWhichCexIsMin) {
               cexRow <- minRowCex
          } else if (rowNumberToUseForSizeCalculation <= rowsAtWhichCexIsMax) {
               cexRow <- maxRowCex
          } else {
               fractionTowardMinimum <- (rowNumberToUseForSizeCalculation - rowsAtWhichCexIsMax)/(rowsAtWhichCexIsMin - rowsAtWhichCexIsMax)
               cexRow <- 0.10 + (maxRowCex-minRowCex)*(1 - (fractionTowardMinimum**0.5))
          }
     }

     if (!is.null(labRow) && numNonBlankRows <= maxNumLabels) {
          ## Right side axis--got to FLIP them around, because of the way we printed the image
          tickLoc.vec <- (0:(nrow(mmm)-1))/(nrow(mmm)-1)
          backwardsLabels.vec <- rev(labRow) ## Gotta reverse the order here!
          ## You can verify that this is correct by examining the values in "heatByCoef" and verifying that
          ## they match up here, too.
          axis(4, at=tickLoc.vec, labels=backwardsLabels.vec, tick=F, las=2, cex.axis=cexRow) # 4 = right axis, usually with gene names
     }

     if (ncol(mmm) == 1) {
          bottomAxisLoc.vec <- c(0.5) ## special case for ncol == 1, so we don't divide by zero!
     } else {
          bottomAxisLoc.vec <- seq(from=0, to=(ncol(mmm)-1)) / (ncol(mmm)-1)
     }
     axis(1, at=bottomAxisLoc.vec, labels=labCol, tick=F, las=2, cex.axis=cexCol) # 1 = bottom axis, usually with array names
     box(lwd=1)

     print.agw("heatmap.agw: Looks like the figure region was an acceptable size to fit the heatmap.")
     # ======================================
     # ==========================================
}








## =============================================================================
## plot.multi.agw: a kind of hack-ish way to plot multiple types of graph from
## the same call. Like if you wanted both a high-res PDF and also a PNG of the same data.
## This is probably not a particularly useful function in most/all cases.

## Example call:

##         plot.multi.agw(type.func=c(
##                        call("pdf", paste(heatmapPath, ".pdf", sep=''), width=15, height=20)
##                        , call("png", paste(heatmapPath, ".png", sep=''), pointsize=18, res=72, width=1800, height=2400))
##                        , plot.func=call("heatmapActual.function")
##                        )
## =============================================================================
plot.multi.agw <- function(type.func, plot.func) {
     ## The arguments to this must be functions generated with "call"
     ## See the examples in the assertions below.
     
     ## type.func needs to be something like "png(filename)" or "pdf(filename)"
     ## plot.func the actual how-to-make-a plot function (example: plot / hist)
     
     for (thePlotMethod in type.func) {
          assert.agw(is.language(thePlotMethod), "The type.func argument must be a VECTOR of types of R functions to output plot data. Example: type.func=c(call(\"png\", width=800), call(\"pdf\", width=10, height=7))")
          assert.agw(is.language(plot.func), "The plot.func argument is the *single* function to call to actually generate the plot. It must be one function, generated with \"call.\" Example: plot.func=call(\"plot\", 1:10, col=\"blue\");")
          print.cyan.agw("plot.multi.agw is running...")
          eval.parent(thePlotMethod) # <-- evaluate in the CALLING environment
          eval.parent(plot.func)     # <-- evaluate in the CALLING environment
          dev.off()
     }
     
}




### ===============================================================================
## Input: a matrix (mm)
## Now, this function goes through each column of the matrix and plots the data from that column as a line.
## If you want to do it by rows instead, transpose the input matrix.
## col can be a VECTOR of colors. So you can pass in something like rainbow( ncol(mm) ) if you want.
## What you end up with is a nice plot with a bunch of lines, like you might see on a sales chart of
## various car brands over time, with several wiggly lines.
## Note that you don't specify the X-axis: it is auto-set to be one integer per Y data point.

## This appears to be basically a reimplementation of "MAPlot," so you should probably
## use that instead.

### ===============================================================================
plotLinesFromMatrixColumnsAgw <- function(mm=NULL, ..., ylim=NULL, col="black", lwd=1, lty=1) {

     ## This is probably a duplication of the built-in function "Matplot"!

     warning("I don't think plotLinesFromMatrixColumnsAgw is probably a good function---you should probably just use the built-in function <matplot> instead!!")
     if (!is.matrix(mm) && is.vector(mm) && (length(mm) > 0)) {
          warning("plotLinesFromMatrixColumnsAgw: the input matrix was actually a one-dimensional vector! Trying to work around this...")
          mm <- matrix( mm , byrow=FALSE, ncol=1)
     }
     assert.agw(!is.null(mm) && is.matrix(mm), "Problem! Input to plotLinesFromMatrixColumnsAgw must be a non-null matrix!")
     if (length(col) < ncol(mm)) { col <- c(col, rep("#00000088", times=ncol(mm))) } ## default fallback color is semi-transparent black
     plot(  x=(1:nrow(mm))
          , y=mm[,1, drop=TRUE]
          , type='l'
          , ylim=ylim
          , ...
          , lwd=lwd, lty=lty, col=col[1])
     for (i in 2:ncol(mm)) {
          lines(  x=(1:nrow(mm))
                , y=(mm[,i, drop=TRUE]), lwd=lwd, lty=lty, col=col[i])
     }
}
### ===============================================================================
### ===============================================================================











### ===============================================================================
### Used in the scatterplots that compare values across two arrays
### ===============================================================================
panel.correlation.local <- function(x, y, digits=2, prefix="", cex.cor=4.0, boxWidth=2) {
     usr <- par("usr"); on.exit(par(usr)) ## restore settings on finishing the plot
     SQUARE_MAX <- 1.0 ; SQUARE_MIN <- 0.0
     MIDDLE_OF_SQUARE <- (SQUARE_MAX + SQUARE_MIN)/2.0
     par(usr = c(SQUARE_MIN, SQUARE_MAX, SQUARE_MIN, SQUARE_MAX))
     #‘usr’ A vector of the form ‘c(x1, x2, y1, y2)’ giving the extremes of the user coordinates of the plotting region.
     
     origR   <- cor(x, y, use="complete", method="pearson")
     absR <- abs(ifelse(is.finite(origR), origR, 0.00))
     MIN_CEX_FAC <- 0.4
     theCex <- cex.cor * max(MIN_CEX_FAC, absR**2) ## don't let the CEX get any smaller than the MIN_CEX

     backgroundColor <- hsv(s = 0.0,  v = 1.0 - 5*(1.0 - max(0.9, absR))) ## uhh... 1.0 = white, and lower scores = a dark-ish gray... kind of hack-ish. Basically it gets as dark as it's going to get around R=0.7 or thereabouts.
     rect(xleft=-100, ybottom=-100, xright=100, ytop=100, col=backgroundColor)
     text(MIDDLE_OF_SQUARE, MIDDLE_OF_SQUARE
          , format(origR, digits=2, nsmall=2, scientific=FALSE)
          , cex=theCex ## <-- The closer R is to zero, the smaller the text size
          )
     box(lwd=boxWidth)
     cat(".") ## progress indicator
}

### ===============================================================================
### Used in the scatterplots that compare values across two arrays
### ===============================================================================
geneST.pairs.plot <- function(..., celKeys) {
     pairsCorMatrixPlot(..., keys=celKeys)
}

### ===============================================================================
## Makes a big multi-plot matrix of correlations between arrays.
## Uses the data in "dataMatrix" to figure out what to plot.
## if you for some reason actually want NO labels, set labelVec to "" in the arguments
### ===============================================================================
pairsCorMatrixPlotAGW <- function(filePath, dataMatrix, labelVec=NULL, keys, main="log2(Intensity).  Red points = within-group comparison. Values in lower left are Pearson's R of the log-transformed values.") {
     
     assert.agw(nrow(keys$"table") == ncol(dataMatrix))
     if (missing(labelVec) || is.null(labelVec)) {
          COLUMN_THAT_PROBABLY_HAS_THE_FILENAMES <- 1 ## usually the first column in the keys file...
          filenamesMostLikely.vec <- keys$table[,COLUMN_THAT_PROBABLY_HAS_THE_FILENAMES]
          labelVec <- paste(filenamesMostLikely.vec, "\n(", keys$groups, ")", sep='')
     }
     labelVec <- gsub("[ ]+", "\n", labelVec, perl=TRUE, ignore.case=TRUE) ## labelVec with a newline added for any spaces
     labelVec <- gsub("\\.(CEL|gpr)", " ", labelVec, perl=TRUE, ignore.case=TRUE )  ## put a newline before the .CEL or .gpr extension at the end of the filename
     labelVec <- gsub("[ _]+vs[.]?[ _]+", "\nvs.\n", labelVec, perl=TRUE, ignore.case=TRUE) ## break up "vs" so it spans multiple lines
     labelVec <- gsub("\n\n+", "\n", labelVec, perl=TRUE, ignore.case=TRUE) ## One newline in a row at most!
     labelVecNoNewlines <- gsub("\n", "", labelVec, perl=TRUE, ignore.case=TRUE)
     
     print.agw("Drawing a \"pairs\" plot matrix with many sub-plots to file <", filePath, ">... (note, this can take a couple of minutes)\n")
     
     #datA <- apply(dataMatrix, APPLY_BY_ROW, mean) ## get the mean of each row
     #datM <- as.matrix(dataMatrix) - datA          ## get the difference of each row

     ## cumul/whichBin: Shows which "bins" each number is in. So like, experiment #2 might be in experimental group #1, if there were 2 replicates of that experiment.
     cumul <- c(0, cumsum(keys$"experimentsInEachGroup")) ## <-- 0 is the leftmost boundary!
     whichBin <- function(value, binMarkers) { ## right sides! The first value is NOT a valid one
          ## Tells you which "bin" an array is in. The bin is the experimental group.
          ## So if you have 9 arrays, in 3 groups (control, control, control, exp1, exp1, exp1, exp2, exp2, exp2),
          ## then the bins are (1, 1, 1, 2, 2, 2, 3, 3, 3)
          ## binMarkers must be sorted in ascending order!
          for (i in 1:(length(binMarkers)-1)) {
               left  <- binMarkers[i]   ;  right <- binMarkers[i+1]
               if (value > left && value <= right) { return(i); }
          }
          return(-1); ## out of bounds
     }

     nGroups <- length(keys$"groups")

     pointAlpha <- 0.75  ## 0.75 = 75% opaque.
     regularPointColor  <- hsv(h=1, s=1, v=0, alpha=pointAlpha) ## The foreground color for points in each "bin"
     backgroundAlpha <- 1.0
     backColors      <- rainbow(n=nGroups*(nGroups-1), s = 0.25, v = 1.0, alpha=backgroundAlpha) ## The background colors for each "bin"
     
     WITHIN_GROUP_POINT_COLOR <- hsv(h=1, s=1, v=1, alpha=pointAlpha)
     WITHIN_GROUP_BACKGROUND_COLOR   <- "white" ## What color are the cells that have the array filenames?

     MINI_PLOT_BORDER_THICKNESS <- 2.0   ## How thick the border around each scatterplot is

     if (ncol(dataMatrix) > 30) {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 100
     } else if (ncol(dataMatrix) > 60) {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 75
     } else {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 200 ## Size of each mini-scatterplot
     }
     MIN_ALLOWED_SIZE_OF_PLOT_IN_PIXELS     <- 1000 ## Minimum size of this plot
     TOTAL_PNG_SQUARE_SIZE <- max((ncol(dataMatrix) * SIZE_IN_PIXELS_PER_SCATTERPLOT), MIN_ALLOWED_SIZE_OF_PLOT_IN_PIXELS)
     # ======================================

     png(filePath, height=TOTAL_PNG_SQUARE_SIZE, width=TOTAL_PNG_SQUARE_SIZE, pointsize=18, res=72)
     par(mar=c("bottom"=2, "left"=2, "right"=2, "top"=8))
     par(cex.axis=1.0, cex.lab=1.2, cex.main=1.0)     
     dimSize <- ncol(dataMatrix)

     whichSubplot <- list("x" = 2, "y" = 1)  ## <-- have to MANUALLY keep track of which sub-plot we are plotting in the  big "upper.panel" function using these vars and the <<- assignment operator

     corrMatrixOutputFilePath <- gsub(".png", ".pearson.corr.matrix.out.txt", filePath)
     print.green.agw("Now writing output correlations to <", corrMatrixOutputFilePath,">...")
     vvv <- matrix(nrow=length(labelVecNoNewlines), ncol=length(labelVecNoNewlines), dimnames=list(labelVecNoNewlines, labelVecNoNewlines))
     for (xx in seq_along(labelVecNoNewlines)) {
          for (yy in seq_along(labelVecNoNewlines)) {
               vvv[xx, yy] <- cor(dataMatrix[,xx], dataMatrix[,yy], use="complete", method="pearson")
          }
     }
     NUM_DECIMAL_PLACES_FOR_CORR_OUTPUT <- 4
     write.table(signif(vvv, NUM_DECIMAL_PLACES_FOR_CORR_OUTPUT), corrMatrixOutputFilePath, sep="\t", col.names=NA, row.names=TRUE, quote=T)
     print.green.agw("[Done] writing output correlations to <", corrMatrixOutputFilePath,">.")
     
     graphics::pairs(dataMatrix
                     , main=main, labels=labelVec
                     , cex.labels=1.5 ## <-- this determines the labels on the diagonal
                     , cex.main=1.0 ## <-- this is related to the size of the graph in a non-intuitive way
                     , lower.panel=function(x,y) {
                          #rect(-10,-10,10,10, col="gray", border=NA) ;
                          cat(paste("Correlating", sep=''))
                          panel.correlation.local(x,y, boxWidth=MINI_PLOT_BORDER_THICKNESS) ;
                          cat("[Done]\n")
                     }
                     , upper.panel=function(x,y) {
                          binX <- whichBin(whichSubplot$x, cumul)
                          binY <- whichBin(whichSubplot$y, cumul)
                          binQuasiRandomIdentifier <- (((binY-1) * nGroups * 13) + binX*17)  ## has to be somehow based on the bins, so that everything in the same bin will be the same color. The exact formula is irrelevant, though, as long as it "seems" random.
                          if (binX == binY) {
                               thePointColor <- WITHIN_GROUP_POINT_COLOR
                               theBackColor <- WITHIN_GROUP_BACKGROUND_COLOR
                          } else {
                               thePointColor <- regularPointColor
                               theBackColor <- backColors[1 + (binQuasiRandomIdentifier %% length(backColors))]
                          }
                          cat(paste("Now calculating the box at X,Y (", whichSubplot$x, ", ", whichSubplot$y, "), which is in bin (", binX, ", ", binY, "), numbered ", binQuasiRandomIdentifier, " and colored <", theBackColor, ">... (", length(x), " points)", sep=''))
                          rect(xleft=min(x)-10, ybottom=min(y)-10, xright=max(x)+10, ytop=max(y)+10, col=theBackColor, border=NA)
                          points(x,y, pch='.', cex=4, col=thePointColor) ## <-- pch='.' is 10 times faster than any other plot character. Use cex=### to adjust the size of the dot.
                          box(lwd=MINI_PLOT_BORDER_THICKNESS)
                          cat(" [Done]\n")
                          whichSubplot$x <<- (whichSubplot$x+1); ## go to the next column... note that this is a GLOBAL assignment! That is important.
                          if (whichSubplot$x > dimSize) {
                               whichSubplot$x <<- (whichSubplot$y+2) ## <-- because upper.panel is a right triangle, we start from a point farther to the right on the plot each time. (This really is supposed to be (xVal <<- yVal+2)!
                               whichSubplot$y <<- (whichSubplot$y+1) ## <-- this only works because the upper.panel is a *right triangle*, so it just happens that the X axis index where we start a new row is the same as the Y axis index. Y axis index is also the same as the X##(whichSubplot$y) #startX;
                          }
                     })
     dev.off()
     print.agw("[Done] with pairs plot.")
}




### ===============================================================================
### 
### ===============================================================================


plot_of_character_frequency_by_position <- function(summarizedCountsFile) {
          # Makes a line plot showing the base frequency by base *position*.
          # Input is normally a "summarized counts" file, which is generated by count_char_freq_at_each_position.pl.
          # This is normally used for RNA-Seq, but can be used for any set of DNA (or other) sequences.

          # summarizedCountsFile is the table that is output by the script "count_char_freq_at_each_position.pl"
          # which is normally run on a fasta file.
          # Note: if you want to get every Nth line from a file, starting with line M:
          #    cat YOURFILE | sed -n 'N~Mp' > OUTPUTFILE

          dat <- read.delim(file=summarizedCountsFile, as.is=TRUE, check.names=FALSE, row.names=1, header=1)

               RSUM_INDEX <- which(rownames(dat) == "RSUM")
               ratios <- apply(dat, 2, function(ccc) { ccc/ccc[RSUM_INDEX] })
               ratios <- ratios[-RSUM_INDEX,]

               LINECOLORS <- rainbow(nrow(ratios))

               if (nrow(ratios) %in% c(4,5) && all(c("A","C","G","T") %in% rownames(dat))) {
                              IS_ACGT_PLOT <- TRUE
                                        theYlab <- "Fraction of each base"
                                        LINECOLORS <- c("#FF0000BB","#000000BB","#0033FFBB","gray","#33CC00BB")
                         } else {
                                        IS_ACGT_PLOT <- FALSE
                                                  theYlab <- "Fraction of each character"
                                   }

               matplot(t(ratios), type='n'
                                    , main=paste("Base frequency at each position in\n\"", summarizedCountsFile, "\"", sep='')
                                    , xlab="Position in each sequence, from the FastQ file", ylab=theYlab
                                    )

               if (IS_ACGT_PLOT) {
                              abline(h=0.25, col="#333333", lty="solid", lwd=2)
                                        abline(h=c(0.2,0.3),  col="#666666", lty="solid" , lwd=1)
                         }
               matplot(t(ratios), type='b', col=LINECOLORS, lty="solid", cex=0.85, pch=rownames(ratios), lwd=2, add=T)
     }



########## ====================== MAIN PROGRAM BELOW ======================

