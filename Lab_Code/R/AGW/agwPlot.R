
## Plotting-related functions.

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
     } else if (grepl("^gr[ea]y", type)) { col <- gray(colRange01) }
     else if (grepl("^brown", type))   { col <- rev(hsv(h=rev(0.2*colRange01), s=colRange01, v=(1.0 - 0.7*colRange01))) }
     else if (grepl("^sepia", type))   { col <- rev(hsv(h=rev(0.3*colRange01), s=colRange01, v=rev(colRange01))) }
     else if (grepl("^heat", type))    { col <- heat.colors(n) }
     else {
          print("Unrecognized type passed into colorGradientAGW:")
          print(type)
          stopifnot(paste("Color type is not recognized. Try something like \"gray.colors\"") == 999)
     }
     if (reverse) { col <- rev(col) }
     return(col)
}




## =================================================================
## Heatmap (by Alex)
## =================================================================
heatmap.agw <- function(m, breaks=12, labRow=rownames(m), labCol=colnames(m), col, main, title="", rowLabelCex=NULL) {
     if (is.vector(m)) {
          m <- as.matrix(m)
     }
     stopifnot(is.matrix(m))
     stopifnot(nrow(m) >= 1)
     stopifnot(ncol(m) >= 1)


     # Generates a three-row multi-part figure.
     # Top part: the caption (title)
     # Middle part: the histogram and distribution key (probably should be made optional, but it's required for now)
     # Bottom part: the heatmap

     m <- t(m)

     stopifnot(breaks >= 3)

     scale01 <- function(x, low = min(x), high = max(x)) { return((x - low)/(high - low)) }

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

     layout(matrix(c(1,2,3), byrow=T, ncol=1), heights=c(lcm(6),lcm(12),1) )

     # ==========================================
     # ======================================

     min.raw <- min(m, na.rm=TRUE)
     max.raw <- max(m, na.rm=TRUE)

     mean.raw <- mean(m, na.rm=TRUE)
     median.raw <- median(m, na.rm=TRUE)


     textDescriptionPlotAGW(paste("Heatmap with ", length(m)
                                  , " values, in ", nrow(m)
                                  , " rows and "
                                  , ncol(m), " columns."
                                  , " Mean = ", format(mean.raw, digits=3, nsmall=1)
                                  , ", median = ", format(median.raw, digits=3, nsmall=1)
                                  , ". ", title
                                  , sep=''), wraplen=120, leftMargin=0.05, topMargin=0.05, cex=1.4)

     # ======================================
     # ==========================================
     # ==========================================
     # ======================================
     par(mar=c("bottom"=3, "left"=9.5, "top"=5, "right"=15.5))

     if (missing(main) || is.null(main)) {
          main <- "Heatmap"
     }
     main <- paste(main, "\n(", nrow(m), " rows by ", ncol(m), " columns)", sep='')

     if (missing(breaks) || is.null(breaks) || (length(breaks) == 1)) {
          ## If breaks was not specified, OR it was a length-one scalar
          breaks <- seq(min.raw, max.raw, length.out=breaks)
     }

     if (missing(col) || is.null(col)) {
          if (totalColors == 2) {
               col <- c("#000066", "#FFFF99")
          } else if (totalColors >= 10) {
               col <- c("#330033", "#440044", "#550055", "#770044", "darkred", heat.colors(totalColors-6))
          } else {
               col <- c("#330033", "#770044", heat.colors(totalColors-2))
          }
          col <- col[1:totalColors]
     } else {
          col <- colors.agw(n = length(breaks)-1, type=col)
     }

     z <- seq(min.raw, max.raw, length = length(col))
     image(z = matrix(z, ncol = 1), col = col, breaks = breaks, xaxt = "n", yaxt = "n")
     par(usr = c(0, 1, 0, 1))
     lv <- pretty(breaks)
     xv <- scale01(as.numeric(lv), min.raw, max.raw)
     axis(1, at = xv, labels = lv)

     h <- hist(m, plot = FALSE, breaks=breaks)
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
     lines(hx, hy/max(hy) * Y_SCALE, lwd=2*HISTOGRAM_LINE_WIDTH, type = "s", col = "white")
     lines(hx, hy/max(hy) * Y_SCALE, lwd=HISTOGRAM_LINE_WIDTH, type = "s", col = "black")

     axis(2, at=pretty(hy)/max(hy) * Y_SCALE, pretty(hy))
     par(cex.main=1.7)
     title(paste("Color Key and Histogram of the Distribution of Heatmap Values\nDashed line = median value (", format(median.raw, digits=3, nsmall=1), ")", sep=''))
     mtext(side=2, paste("Count (out of ", length(m) ,")", sep=''), line=3)
     box(lwd=1)
     # ======================================
     # ==========================================
     # ==========================================
     # ======================================
     par(mar=c("bottom"=25, "left"=8, "top"=8, "right"=20), cex.main=2.2)
     par(cex = 0.5) ## Expansion factor
     image(m, axes=F, main=main, col=col)

     numRows <- ncol(m) ## -- yes, it really is NCOL, because we transposed m due to the way "image" draws things
     numCols <- nrow(m) ## -- yes, it really is NROW, because we transposed m due to the way "image" draws things
     if (is.null(rowLabelCex)) {
          rowLabelCex = 1.0

          rowsAtWhichCexIsMin <- 800 ## The number of rows at which the text is the tiniest (or more rows than this)
          rowsAtWhichCexIsMax <- 100 ## The number of rows at which the text is the biggest (or fewer rows than this)
          stopifnot(rowsAtWhichCexIsMin > rowsAtWhichCexIsMax)
          maxRowCex <- 1.00
          minRowCex <- 0.10

          if (numRows >= rowsAtWhichCexIsMin) {
               rowLabelCex <- minRowCex
          } else if (numRows <= rowsAtWhichCexIsMax) {
               rowLabelCex <- maxRowCex
          } else {
               fractionTowardMinimum <- (numRows - rowsAtWhichCexIsMax)/(rowsAtWhichCexIsMin - rowsAtWhichCexIsMax)
               rowLabelCex <- 0.10 + (maxRowCex-minRowCex)*(1 - (fractionTowardMinimum**0.5))
          }
     }

     MAX_NUM_ROWS_TO_LABEL <- 1000
     if (numRows <= MAX_NUM_ROWS_TO_LABEL) {
          axis(4, at=(0:(ncol(m)-1))/(ncol(m)-1), labels=labCol, tick=F, las=2, cex.axis=rowLabelCex) # 4 = right axis, usually with gene names
     }

     axis(1, at=(0:(nrow(m)-1))/(nrow(m)-1), labels=labRow, tick=F, las=2, cex.axis=1.5) # 1 = bottom axis, usually with array names
     box(lwd=1)
     # ======================================
     # ==========================================
}
