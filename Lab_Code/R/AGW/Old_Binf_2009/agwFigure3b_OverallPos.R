#
#
##    source("/Users/alexgw/R/src/agwFigure3b_OverallPos.R");   ## <-- To reload this file

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

# Generate the aggregate plots showing positive control recovery
# rate for ALL modules.

DIR_3B = "Figure3b_OverallPosRecovery"


FINALIZE = TRUE # <-- if this is set, makes a stripped-down "final" figure. Otherwise plots a more-informative-but-also-confusing full figure
REQUIRED_NUM_HITTING_DRUGS <- 4 # require at least this many hitting drugs for a module before we include it

SHOULD_OVERPLOT_LABELS = TRUE #!(FINALIZE) # whether or not we label each pathway


agwSubAggregateLinePlot <- function(dataCollection, xx, yy, nn=NULL, xlim=NULL, ylim=NULL, ...) {
     # Plot only one single line plot. Note: requires the PNG stuff to be set up previously!
     
     #browser()
     stopifnot(length(xx) == length(yy)) ## <-- X and Y vector lengths to aggregateLinePlots must match
     
     finalYLimits = NA
     if (!SHOULD_OVERPLOT_LABELS) { finalYLimits = c(-1.0, 1.0) ; }
     else { finalYLimits = c(-1.2, 1.2) }
     if (!is.null(ylim)) {  finalYLimits  <-  ylim  }
     
     if ((length(xx) == 0) || (length(yy) == 0)) {
          stop("ERROR---no valid data to agwSubAggregateLinePlot");
     }
     
     finalXLimits = NA
     if (FINALIZE) {
          finalXLimits = c(0, agwMax(sapply(xx, function(a) { max(a, na.rm=TRUE) } )))
     } else {
          finalXLimits = c(0, agwMax(sapply(xx, function(a) { max(a, na.rm=TRUE) + 5 } )))
     }
     if (!is.null(xlim)) {  finalXLimits  <- xlim  }
     
     lineColors    <-  agwColorsFromStrings(paste(names(xx), "tmp123"))
     lineNames     <-  agwAbridgeModuleNames(names(xx))
     
     plotItemScale <-  1.0 #25
     
     zeroLineColor      <- "dark gray"
     hlineColor         <- "gray"
     vlineColor         <- "gray"
     
     left1Color         <- "#00000055"
     left2Color         <- "#00000022"
     outOfBoundsOpaque  <- "#999999FF"
     worseThanZeroColor <- "#44444455"
     
     zeroLineWidth      <-  3
     markLineWidth      <-  3
     smoothedLineWidth  <-  2
     
     plot(x=c(1), y=c(1), type='n', xlab="PGI Score", ylab="Recovery rate of ORF differential sensitivity"
          , xlim=finalXLimits, ylim=finalYLimits, ...) ## <- start the plot...
     rect(xleft=0, ybottom=-10, xright=100, ytop=0, col=worseThanZeroColor, border=0) ## show the "worse than 0" region

     abline(v=seq(0, by=5, length.out=20), col=vlineColor, lwd=markLineWidth)
     
     ##abline(h=0, col=zeroLineColor, lwd=zeroLineWidth*plotItemScale)
     #rect(xleft=0, ybottom=-10, xright=2,  ytop=10, col=left1Color, border=0)
     #rect(xleft=2, ybottom=-10, xright=5,  ytop=10, col=left2Color, border=0)
     #for (i in 1:20) {
     #     rect(xleft=(i*10), ybottom=-10, xright=(i*10)+5,  ytop=10, col=left2Color, border=0)
     #}
     
     ## Mark the "out of bounds" region above/below +/- 1
     rect(xleft=0, ybottom=-10, xright=100,  ytop=-1 , col=outOfBoundsOpaque, border=0)
     rect(xleft=0, ybottom=1,   xright=100,  ytop=100, col=outOfBoundsOpaque, border=0)
     ##abline(h=c(-1,1), col="black", lwd=zeroLineWidth*plotItemScale)
     
     abline(v=0, col="black", lwd=zeroLineWidth*plotItemScale)
     
     ## The horizontal lines
     abline(h=c(-0.75,-0.5,-0.25,0.25,0.5,0.75), col=hlineColor, lwd=markLineWidth*plotItemScale, lty="dotted")
     
     for (i in 1:length(xx)) {
          modName <- (names(xx))[i]
          numDrugsThatTargetThisModule <- length(dataCollection[[modName]])
          #print(numDrugsThatTargetThisModule)
          
          points(    x = xx[[i]]
                 ,    y = yy[[i]]
                 ,  col = lineColors[i]
                 ,  lwd = (smoothedLineWidth*plotItemScale)
                 , type = 'o'
                 ,  pch = PCH.SMALL.CIRCLE
                 ,  cex = (0.95 * plotItemScale)
                 )
          
          if (!FINALIZE) {
               if (!is.null(nn) && length(nn[[i]]) > 0) {
                    thresholdForShowingNumberOfDataPoints <- 9 ## <-- if this many (or fewer) points were used to make a smoothed line, show it!
                    pchAboveThresh <- ''
                    text(xx[[i]], yy[[i]]
                         , labels = ifelse(nn[[i]] <= thresholdForShowingNumberOfDataPoints, nn[[i]], pchAboveThresh)) ## Plot a text indicator of the number of data points that were used to compute this average
               }
          }
          
     }
     
     if (SHOULD_OVERPLOT_LABELS) {
          for (i in 1:length(xx)) {
               modName <- (names(xx))[i]
               numDrugsThatTargetThisModule <- length(gvDrugTargetMapping[[modName]])
               
               singleText <- function(whichIndex) {
                    ## Draw the module name associated with the line.
                    ## The index is "which point on the line should we use for labeling the entire line?"
                    ## Typically this is the first or last point, but a middle point is also OK.
                    yloc = yy[[i]][whichIndex]
                    baseRotationAngle = 20
                    rotateAmt = baseRotationAngle
                    if (!is.na(yloc) && length(yloc) > 0 && yloc < 0.25) {
                         rotateAmt = -baseRotationAngle
                    }
                    par(srt = rotateAmt) ## Rotate the text ("String rotate")
                    if (length(xx[[i]][whichIndex]) > 0) {
                         theLabel <- agwPerlSub("KNOWN:DB_FROM_HUMAN:","DBANK:"
                                                , agwPerlSub("(KNOWN:HAND_CURATED:)","HAND:",lineNames[i]))
                         ## Label the drug lines
                         text(  xx[[i]][whichIndex]
                              , yloc
                              , labels=agwGlue("<---------------- ", theLabel
                                , " (", numDrugsThatTargetThisModule, " hitters, ", round(gvCompleteCollection[[i]]$moduleBaselineFract,2), ", ", round(gvCompleteCollection[[i]]$moduleBaselineFract,2), " baseline)")
                              , cex=(0.9 * plotItemScale)
                              , pos=POS.RIGHT
                              )
                    }
               }
               
               singleText(which.max(yy[[i]]))
               singleText(which.max(xx[[i]]))
               ##singleText( length(xx[[i]]) )
               ##singleText( ceiling(length(xx[[i]]) / 2) )
               ##singleText( 1 )
          }
     }
     
     par(srt = 0) ## Rotate the text back to the default horizontal ("String rotate")
     legend("bottomright"
            , legend=lineNames
            , fill=lineColors
            , cex=(0.75*plotItemScale)
            , bty='n'
            )
}

agwPlotAggregateLinePlots <- function(dataCollection, filetype, directory, file, finalize, xx, yy, yybase=NULL, nn=NULL, xlim=NULL, ylim=NULL, main1=NULL, main2=NULL, ...) {
     # Plot two line plots--one showing the absolute results, one showing the normalized results.
     ## xx and yy are the "here are the points that were NOT adjusted"
     ## yybase is "the yy points MINUS some kind of baseline"
     agwPreparePlot(  directory = directory
                    , filetype  = filetype
                    ,      file = file
                    ,    width  = 4*STANDARD.PLOT.WIDTH
                    ,    height = 5*STANDARD.PLOT.HEIGHT
                    ,       res = 140
                    , pointsize = (DEFAULT.POINTSIZE*1.2))
     
     stackSplitLoc = 0.5
     split.screen(matrix(data=c(
                         0.0, 1.0, stackSplitLoc, 1.00, ## left, right, bottom, top
                         0.0, 1.0, 0.00, stackSplitLoc), byrow=T, ncol=4))
     
     screen(1) ;  agwSubAggregateLinePlot(dataCollection=dataCollection, xx=xx, yy=yy, nn=nn, xlim=xlim, ylim=ylim, main=main1)
     if (!is.null(yybase)) {
          screen(2) ;  agwSubAggregateLinePlot(dataCollection=dataCollection, xx=xx, yy=yybase, nn=nn, xlim=xlim, ylim=ylim, main=main2)
     }
     close.screen(all=TRUE)
     agwFinishPlot()
}


generateAggregatePlots <- function(ftype, shouldLogXAxis, shouldMakeFinalForPaper=FALSE) {
     agwGlobalLoad("gvCompleteCollection") # <-- set in agwFigureA_DetailPlots.R, although the immediate call is in agwDrugCalcFunctions.R
     enough <- sapply(gvDrugTargetMapping
                      , function(a) { (length(a) >= REQUIRED_NUM_HITTING_DRUGS) } # are there enough hitting drugs for us to be happy with this module?
                      )
     modNamesInCompleteCollectionWithEnoughHittingDrugs  <- intersect( names(gvCompleteCollection) , names(gvDrugTargetMapping[enough]))
     hasEnoughHittingDrugsCollection <- gvCompleteCollection[modNamesInCompleteCollectionWithEnoughHittingDrugs]  #[enough] ## only the modules with a sufficient number of drugs (as determined by REQUIRED_NUM_HITTING_DRUGS)
     #browser()

     recoveryPrefix = "Recovery"
     if (shouldLogXAxis) { recoveryPrefix = agwGlue("Log10_X_Recovery") }
     
     processXValues <- function(xVec) {
          if (shouldLogXAxis) {
               return(log(xVec+1, base=2)) # 0 used to be the worst, so now we make it 1
          } else {
               return(xVec)
          }
          
     }
     
     agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                               , filetype = ftype
                               , directory = agwGlobalOutPath(DIR_3B)
                               , finalize = shouldMakeFinalForPaper
                               , file = agwGlue(recoveryPrefix)
                               ,     xx=lapply(hasEnoughHittingDrugsCollection, function(a) { processXValues(a$smoothedWide$x) })
                               ,     yy=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$y })
                               ,     nn=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$n })
                               , yybase=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$y - a$moduleBaselineFract })
                               ,  main1=agwGlue("Raw recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ". Only modules with at least ", REQUIRED_NUM_HITTING_DRUGS, " hitting drugs.")
                               ,  main2=agwGlue("Better-than-baseline Recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus SL baseline (i.e., the rate obtained by guessing an in-module gene.")
                               )

     if (!shouldLogXAxis) {
          # the log-scale-on-x-axis plots do not need a "zoomed" version
          agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                                    , filetype = ftype
                                    , directory = agwGlobalOutPath(DIR_3B)
                                    , finalize = shouldMakeFinalForPaper
                                    , file = agwGlue(recoveryPrefix, "_ZOOM")
                                    ,     xx=lapply(hasEnoughHittingDrugsCollection, function(a) { processXValues(a$smoothedWide$x) })
                                    ,     yy=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$y })
                                    ,     nn=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$n })
                                    , yybase=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedWide$y - a$moduleBaselineFract })
                                    ,  main1=agwGlue("Smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ".")
                                    ,  main2=agwGlue("Smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus module baseline (i.e., the rate obtained by guessing an in-module gene).")
                                    ,   xlim=c(0,7)
                                    )
     }
     
     agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                               , filetype = ftype
                               , directory = agwGlobalOutPath(DIR_3B)
                               , finalize = shouldMakeFinalForPaper
                               , file = agwGlue(recoveryPrefix, "_Narrow_Smoothing")
                               ,     xx=lapply(hasEnoughHittingDrugsCollection, function(a) { processXValues(a$smoothedNarrow$x) }) # rawUnsmoothed
                               ,     yy=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$y })
                               ,     nn=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$n })
                               , yybase=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$y - a$moduleBaselineFract })
                               ,  main1=agwGlue("Narrow smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ".")
                               ,  main2=agwGlue("Narrow smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus module baseline (i.e., the rate obtained by guessing an in-module gene).")
                               )

     if (!shouldLogXAxis) {
          # the log-scale-on-x-axis plots do not need a "zoomed" version
          agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                                    , filetype = ftype
                                    , directory = agwGlobalOutPath(DIR_3B)
                                    , finalize = shouldMakeFinalForPaper
                                    , file = agwGlue(recoveryPrefix, "_Narrow_Smoothing_ZOOM")
                                    ,     xx=lapply(hasEnoughHittingDrugsCollection, function(a) { processXValues(a$smoothedNarrow$x) }) # rawUnsmoothed
                                    ,     yy=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$y })
                                    ,     nn=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$n })
                                    , yybase=lapply(hasEnoughHittingDrugsCollection, function(a) { a$smoothedNarrow$y - a$moduleBaselineFract })
                                    ,  main1=agwGlue("Narrow smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ".")
                                    ,  main2=agwGlue("Narrow smoothed recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus module baseline (i.e., the rate obtained by guessing an in-module gene).")
                                    
                                    ,   xlim=c(0,7)
                                    )
     }
     
     cxx <- list() # <-- these are the cumulativeRecov data collections, but with all of the (x = 0) (i.e., pgi score of 0)
     cyy <- list() #     ORFs collapsed into one single data point instead of like.... 4000 or so.
     cnn <- list()
     cyyRelative <- list()
     
     for (i in 1:length(hasEnoughHittingDrugsCollection)) {
          # Here we go through and merge all the data points at x = 0.
          # Because there are SO many ORFs that have no PGI/SGI/synthetic lethal overlap with 
          x1 <- hasEnoughHittingDrugsCollection[[i]]$cumulativeRecov$x
          y1 <- hasEnoughHittingDrugsCollection[[i]]$cumulativeRecov$y
          n1 <- hasEnoughHittingDrugsCollection[[i]]$cumulativeRecov$n
          
          cxx[[i]] <- x1
          cyy[[i]] <- y1
          cnn[[i]] <- n1
          
          whichZeroX <- which(x1 == 0)
          if (length(whichZeroX) > 0) {
               firstZeroLoc <- min(whichZeroX) # the first value that is zero. The input array better be sorted, though!
               newYmean <- mean(y1[x1 == 0])     # All the SUBSEQUENT values in x1 are expected to be zero after firstZeroLoc
               newNmax  <- max(n1)
               if (firstZeroLoc >= 3) {
                    cxx[[i]] <- c(x1[1:(firstZeroLoc-1)], 0)
                    cyy[[i]] <- c(y1[1:(firstZeroLoc-1)], newYmean)
                    cnn[[i]] <- c(n1[1:(firstZeroLoc-1)], newNmax)
               } else {
                    # don't set newX, etc...
               }
               
          } else {
               # don't set newX, etc...
          }
          
          cyyRelative[[i]] <- cyy[[i]] - hasEnoughHittingDrugsCollection[[i]]$moduleBaselineFract  # <-- y points minus base
     }
     names(cxx) = names(hasEnoughHittingDrugsCollection)
     names(cyy) = names(hasEnoughHittingDrugsCollection)
     names(cnn) = names(hasEnoughHittingDrugsCollection)
     names(cyyRelative) = names(hasEnoughHittingDrugsCollection)
     
     #browser()
     
     agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                               , filetype  = ftype
                               , directory = agwGlobalOutPath(DIR_3B)
                               , finalize  = shouldMakeFinalForPaper
                               , file = "Recovery_Cumulative"
                               ,     xx=cxx
                               ,     yy=cyy
                               ,     nn=cnn
                               , yybase=cyyRelative
                               ,  main1=agwGlue("Raw cumulative recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ".")
                               ,  main2=agwGlue("Cumulative recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus module baseline (i.e., the rate obtained by guessing an in-module gene).")
                               )
     
     agwPlotAggregateLinePlots(dataCollection = hasEnoughHittingDrugsCollection
                               , filetype  = ftype
                               , directory = agwGlobalOutPath(DIR_3B)
                               , finalize  = shouldMakeFinalForPaper
                               , file = "Recovery_Cumulative_ZOOM"
                               ,     xx=cxx
                               ,     yy=cyy
                               ,     nn=cnn
                               , yybase=cyyRelative
                               ,  main1=agwGlue("Raw cumulative recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ".")
                               ,  main2=agwGlue("Cumulative recovery rates of diff. sensitivity, at sens ratio threshold of ", foldDifferenceForSensitivityToCount, ", minus module baseline (i.e., the rate obtained by guessing an in-module gene).")
                               ,   xlim=c(0,7)
                               )
     
}


agwFigure3b_OverallPos_main <- function() {
     #generateAggregatePlots(ftype="png");
     generateAggregatePlots(ftype="pdf", shouldLogXAxis=FALSE, shouldMakeFinalForPaper=FALSE);
     generateAggregatePlots(ftype="pdf", shouldLogXAxis=TRUE,  shouldMakeFinalForPaper=FALSE);
}



##
