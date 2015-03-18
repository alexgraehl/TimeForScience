#
#
## source("/Users/alexgw/R/src/agwFigureA_DetailPlots.R");     ## <-- To reload this file


source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");


agwFigureA_GenerateManyDetailPlots <- function() {
     ##par(ask=FALSE)
     ##par(ask=TRUE)
     agwGlobalLoad("gvDrugTargetMapping")
     agwGlobalLoad("gvNeighborStatsList")
     agwGlobalLoad("gvMemberStatsList")

     for(i in 1:length(gvDrugTargetMapping)) {
          name <- names(gvDrugTargetMapping)[i]
          agwSensVsOverlapPlot(moduleName=name
                               , theStatCollection=gvNeighborStatsList
                               , theModuleMemberCollection=gvMemberStatsList
                               , filePath=agwGlobalOutPath("FigureA_Indicators/")
                               , fileName=agwGlue("Indicators_for_",name))
          ##readline() ## <-- pauses after making ANY plot (good for interactive plot viewing)
     }
     ##par(ask=FALSE)
     
     if (!agwHasContent("gvCompleteCollection") || length(gvCompleteCollection) == 0) {
     	stop("Huh! For some reason, gvCompleteCollection failed to get ANY elements added to it. Probably some paths are wrong somewhere!")	
     }
}



agwSensVsOverlapPlot <- function(moduleName, theStatCollection, theModuleMemberCollection, filePath=NULL, fileName=NULL) {
     ## "gvCompleteCollection" is set here
     slList     <- theStatCollection[[moduleName]]
     memberList <- theModuleMemberCollection[[moduleName]]
     theModTargetsThisManyORFs <- length(slList)
     agwGlobalLoad("gvPathwayMembershipHash")
     
     orfsInThisModuleList <- agwHashGet(hash=gvPathwayMembershipHash, key=moduleName)
     ## (Higher in the original scale = worse growth)
     ## TARGETS it - (doesn't target it)


     uninterestingCutoff <- 2.0
     interestingCutoff   <- 5.0
     uninterestingOverlapColor       <- "#DDEE99"
     slightlyInterestingOverlapColor <- "#EEFFCC"

     densityFillColor     <- "#B0FFA0"
     darkDensityFillColor <- "#229922"
     ptOutlineColor       <- "gray"
     ptFillColor          <- "black"

     linearFitColor       <- "brown"
     linearFitLwd         <- 3

     thresholdBarColor    <- "red"  ## the "ORFs must be *this* differentially sensitive" threshold
     thresholdBarLwd      <- 2

     baselineModuleColor  <- "#00114490" ## <-- should be somewhat transparent, as it over-plots!
     baselineSLColor      <- "#00ee00AA"
     baselineLwd          <- 3
     baselineLty          <- "solid"

     yesThreshColor       <- "#EE3333EE"
     noThreshColor        <- "#0000FF99" ## blue-ish, semi-transparent

     smallLabelCex      <- 0.6
     pointScalingFactor <- 0.8 ## <-- size of the triangles

     highScoringCex = pointScalingFactor * 1.5 ## the size of the top-scoring gene points
     highScoringLwd = 3  ## the line width for the top-scoring genes (by SL overlap with the target pathway)

     inModulePtFillColor   <- "#00FF00"
     inModulePtBorderColor <- "#006600"

     inModulePtScale       <- pointScalingFactor * 2 ## the size of the points for the in-module ORFs
     inModulePtLWD         <- 2  ## the width of the border for the points of ORFs that are *in* the drug-targeted module

     smoothedLineWidth     <- 3  ## <-- the widths for the "smoothed" lines on the bottom plot
     ##agwSmoothNarrowColor  <- "orange"
     agwSmoothWideColor    <- "red"
     agwSmoothAccumColor   <- "black"
     agwAccumBaseColor     <- "purple"
     superSmoothedColor    <- "blue"
     loeSmoothedColor      <- "green"

     yesSuper <- FALSE ## <-- show the "super smoothed" curve?
     yesLoess <- FALSE ## <-- show the loess-smoothed curve?

     geneLabelDisplayThresholdY <- 1 ## <-- display any gene labels with abs(y) values that are > this amount
     geneLabelDisplayThresholdX <- uninterestingCutoff ## <-- display any gene labels with x values that are > this amount

     legendCex = 1.0

     linePtCex = 0.55 ## <-- size of the "nubs" on the lines

     zeroLineLwd = 1
     thinDelimiterLineLwd = 1
     thinDelimiterLineStyle = "dotted"

     topNMaxVec    <- c(5,10)  ## How many "top results" are we fishing for? This is a VECTOR!. If you add another element here, you get more "top overlap" values being displayed as highlighted points
     highScoringColorsVec <- c("purple","yellow")

     slDiffVec <- as.vector(calcManyDiffSens(slList), mode="double")
     slOverVec <- as.vector(sapply(slList, "[[", "overlap"), mode="double") # it's -log10(pval) (0 is worst, 1000 would be ridiculously good)

     memberDiffVec <- calcManyDiffSens(memberList)
     memberOverVec <- as.vector(sapply(memberList, "[[", "overlap"), mode="double")

     agwSmoothWindowSize <- 0.50
     agwSmoothNumSteps   <- 50

     ##   print("---")
     ##   agwPrint("Member list: ", memberList)
     ##   print(length(memberOverVec)) ; print(memberOverVec)
     ##   print(length(memberDiffVec)) ; print(memberDiffVec)
     ##   print("---")

     isSensBooleanVec <- isDifferentiallySensitive(slDiffVec)
     isSensColorVec   <- ifelse(isSensBooleanVec == TRUE, yesThreshColor, noThreshColor)
     isSensShapeVec   <- ifelse(isSensBooleanVec == TRUE, PCH.BORDER.UP.TRIANGLE, PCH.BORDER.DOWN.TRIANGLE)
     isSensVec        <- ifelse(isSensBooleanVec == TRUE, 1, 0)
     rm(isSensBooleanVec)

     minLength <- agwMin(length(na.omit(slDiffVec)), length(na.omit(slOverVec)))
     requiredNumberOfNonNADataPoints <- 1 ## Spearman needs 2, Pearson requires 3, density requires... 3 I think
     if (minLength < requiredNumberOfNonNADataPoints) {
          fileName <- agwGlue("Z_Empty_",fileName)
     }

     LARGE.MULTI.PLOT.WIDTH  = 3*STANDARD.PLOT.WIDTH
     LARGE.MULTI.PLOT.HEIGHT = 2*STANDARD.PLOT.HEIGHT
     
     agwPreparePlot(directory = filePath
                    , file = fileName
                    , width = LARGE.MULTI.PLOT.WIDTH
                    , height = LARGE.MULTI.PLOT.HEIGHT
                    , res = HIGH.RES.DPI
                    , pointsize = DEFAULT.POINTSIZE)

     if (minLength < requiredNumberOfNonNADataPoints) {
          ## Plots a big "sorry, this won't work" graph
          notEnoughStr <- agwGlue("Not enough data to plot\n<",moduleName,">, (",length(slList), ", ", length(na.omit(slDiffVec))," and ", length(na.omit(slOverVec)), ")")
          plot(c(1),c(1),pch="X",cex=10,col="dark gray",main=notEnoughStr,ylab="",xlab=notEnoughStr,axes=FALSE)
          print(notEnoughStr)
          agwFinishPlot()
          return() ## Not enough data to generate a meaningful plot
     }

	 if (!agwHasContent("gvCompleteCollection")) {
          gvCompleteCollection <<- list() ## <-- A collection of the "smoothed" lines showing differentially sensitivity vs. hypergeometric overlap
     }
     
     numSensitiveOverall <- length(na.omit(slDiffVec[isDifferentiallySensitive(slDiffVec)]))
     numValidDiffPoints  <- length(na.omit(slDiffVec))

     numSensitiveInModule <- length(na.omit(memberDiffVec[isDifferentiallySensitive(memberDiffVec)]))
     numValidInModule   <- length(na.omit(memberDiffVec))

     topIndicesByOverlapVec <- order(sort(slOverVec,decreasing=FALSE))
     topByOverlapList       <- slList[topIndicesByOverlapVec ]
     numCorrectVec <- numSkippedVec <- numCheckedVec <- rep(0,length(topNMaxVec))

     stopifnot(length(topNMaxVec) == length(highScoringColorsVec))

     for (n in 1:length(topNMaxVec)) {
          checkMaxNum <- agwMin(topNMaxVec[n],length(na.omit(slDiffVec)),length(na.omit(slOverVec)))
          for (i in 1:length(slList)) {
               if (numCheckedVec[n] >= checkMaxNum) { break; }
               item <- slList[[i]]
               if (is.na(slDiffVec[i])) {
                    numSkippedVec[n] <- (numSkippedVec[n]+1) ## Skip points without experimental data!
                    next; ## next iteration of the loop...
               }

               #agwPrint("ORF ",item$orf," (",item$gene,"), ranked #",i," in SGA overlap (-log10(p)=",round(item$overlap,1),") with the module <",moduleName,"> has differential sensitivity ",round(slDiffVec[i],2))
               if (isDifferentiallySensitive(slDiffVec[i])) {
                    numCorrectVec[n] <- (numCorrectVec[n]+1)
                    #agwPrint("So it is differentially sensitive!")
               }
               numCheckedVec[n] <- (numCheckedVec[n]+1)
          }
          overallAccuracyOfTop <- (numCorrectVec[n]/numCheckedVec[n])
          #agwPrint("Overall accuracy looking at just the top ",checkMaxNum," out of ",numCheckedVec[n]," tested ORFs: ", numCorrectVec[n],"/",numCheckedVec[n]," (",round(overallAccuracyOfTop,2),") (Skipped ", numSkippedVec[n]," ORFs with high overlap but no experimental data out of the entire pool of ",length(slList),")")
     }


     print(agwGlue("Figure A: Working on module ", moduleName))
     print(length(na.omit(slOverVec))) ;  print(length(na.omit(slDiffVec)))

     hasValueInBothVec <- !is.na(slDiffVec) & !is.na(slOverVec)
     hasValidDataCount <- length(slDiffVec[ hasValueInBothVec ]) ## Number of things that are not NA
     stopifnot(length(slDiffVec) == length(slOverVec))

     inModulePercentRight <- 100*numSensitiveInModule/numValidInModule

     baselinePercentRight <- 100*numSensitiveOverall/numValidDiffPoints
     recoveredPercentVec  <- 100*numCorrectVec/numCheckedVec

     numDrugsThatTargetModule <- 0
     if (length(slList) >= 1) {
          numDrugsThatTargetModule <- slList[[1]]$nTargetingModule
          numDrugsNotTargetingModule <- slList[[1]]$nNotTargetingModule
     }

     pearCor <- NULL ; if (minLength >= 3) {
          suppressWarnings(pearCor <- cor.test(slOverVec, slDiffVec, method="pearson"))
     }
     spearCor <- NULL ; if (minLength >= 2) {
          suppressWarnings(spearCor <- cor.test(slOverVec, slDiffVec, method="spearman"))
     }
     linearFit <- NULL ; if (minLength >= 3) {
          suppressWarnings(linearFit <- lm(slDiffVec ~ slOverVec))
     } ## Line of best fit / linear fit

     stackSplitLoc = 0.4
     split.screen(matrix(data=c(
                         0.0, 1.0, stackSplitLoc, 1.00, ## left, right, bottom, top
                         0.0, 1.0, 0.00, stackSplitLoc), byrow=T, ncol=4))

     pearsonPrint  <- NULL ; if (!is.null(pearCor)) { pearsonPrint <- round(pearCor$estimate,2) }
     spearmanPrint <- NULL ; if (!is.null(spearCor)) { spearmanPrint <-round(spearCor$estimate,2) }

     xHypergeoLabel="-log10(p) Hypergeometric score for the ORF's SGI overlap with the module in question"
     plotYlim <- mirrorYlim(slDiffVec)
     plotXlim <- c(0, agwMax(5, agwMax(slOverVec)+1))
     screen(1) ; plot(x=slOverVec, y=slDiffVec
                      , ylim=plotYlim
                      , xlim=plotXlim
                      , xlab=xHypergeoLabel
                      , ylab="Differential Sensitivity: (med(hitting)-med(other))/(median.abs.dev)"
                      , main=agwGlue("About <",agwAbridgeModuleNames(moduleName), "> (a module of size ", length(orfsInThisModuleList), ")", "\n",theModTargetsThisManyORFs," ORFs with nonzero SGI overlap (",hasValidDataCount
                        , " with valid data). Pearson=", pearsonPrint
                        , ", Spearman=", spearmanPrint
                        )
                      , pch=' ' ## <-- don't ACTUALLY plot anything visible! This just sets up the window
                      )
     
     rect(xleft=0, ybottom=agwMin(plotYlim), xright=uninterestingCutoff+0.1, ytop=agwMax(plotYlim), col=uninterestingOverlapColor, border=0) ## Plot the "uninteresting" areas
     rect(xleft=uninterestingCutoff, ybottom=agwMin(plotYlim), xright=interestingCutoff, ytop=agwMax(plotYlim), col=slightlyInterestingOverlapColor, border=0)
     
     if (minLength >= 3) {
          agwPlotLinesViolinDensity(slOverVec, col=densityFillColor, scaleToY=plotYlim[2], border=darkDensityFillColor)
     }
     
     if (!is.null(linearFit) && !is.null(linearFit$a) && !is.null(linearFit$b)) {
          abline(linearFit, col=linearFitColor, lwd=linearFitLwd, lty="dashed") ## line of best fit
     }
     abline(h=0,col="black", lwd=zeroLineLwd)
     
     ## Show the differential-sensitivity threshold
     abline(h=sensThresh, col=yesThreshColor, lwd=thresholdBarLwd)
     text(x=plotXlim[2],y=sensThresh+0.12,labels=agwGlue("Significant differential sensitivity threshold (",round(sensThresh,2),")"),col=thresholdBarColor,cex=smallLabelCex,pos=POS.LEFT)
     
     ## Draw the gene names for the synthetic lethal overlap-chosen points
     if (length(slList) > 0) {
          validsVec <- (slOverVec > geneLabelDisplayThresholdX) | (abs(slDiffVec) > geneLabelDisplayThresholdY)
          if (sum(validsVec) > 0) {
               text(      x = slOverVec[validsVec]
                    ,     y = slDiffVec[validsVec],
                    , labels=sapply(slList, "[[", "gene")[validsVec ]
                    , pos=POS.RIGHT, col="black", cex=smallLabelCex) ## Gene names
          }
          points(x=slOverVec, y=slDiffVec, cex=pointScalingFactor
                 , col=isSensColorVec, bg=isSensColorVec, pch=isSensShapeVec)
     }
     
     if (length(memberList) > 0) {
          validsVec <- (memberOverVec > geneLabelDisplayThresholdX) | (abs(memberDiffVec) > geneLabelDisplayThresholdY)
          if (sum(validsVec) > 0) {
               ## Draw the gene names for the in-module points
               text(       x = memberOverVec[validsVec]
                    ,      y = memberDiffVec[validsVec]
                    , labels = sapply(memberList,"[[","gene")[validsVec]
                    ,    pos = POS.RIGHT, col=inModulePtBorderColor, cex=smallLabelCex) ## Gene names
          }
     }

     ## Specially-plot the highest scoring (by overlap score) genes
     for (n in length(topNMaxVec):1) {
          numTopXPointsToPlot <- topNMaxVec[n]
          points(  x=slOverVec[1:numTopXPointsToPlot]
                 , y=slDiffVec[1:numTopXPointsToPlot]
                 , pch=isSensShapeVec, col="black", bg=highScoringColorsVec[n]
                 , cex=highScoringCex , lwd=highScoringLwd )
     }

     ## Draw the in-module points
     if (length(memberList) > 0) {
          points(x=memberOverVec, y=memberDiffVec, cex=inModulePtScale
                 , col=inModulePtBorderColor, bg=inModulePtFillColor, pch=PCH.DIAMOND, lwd=inModulePtLWD)
     }

     legend("topright"
            , legend=c(agwGlue("No. ORFs with sensitivity data: ", hasValidDataCount)
              , agwGlue("No. drugs/replicates targeting this module: ", numDrugsThatTargetModule)
              , agwGlue("No. drugs/replicates targeting other modules: ", numDrugsNotTargetingModule)
              , agwGlue("Total drugs/replicates with known targets: ", (numDrugsThatTargetModule+numDrugsNotTargetingModule))
              , agwGlue("Pearson: ",pearsonPrint,", Spearman: ",spearmanPrint)
              )
            , cex=legendCex, box.lwd=0, bty='n')

     legend("bottomleft"
            , legend=c("Density (scaled)"
              , agwGlue("Overlap score < ", uninterestingCutoff)
              , agwGlue("Overlap score < ", interestingCutoff)
              , agwGlue("Line of best fit"))
            , fill=c(densityFillColor, uninterestingOverlapColor, slightlyInterestingOverlapColor, linearFitColor)
            , cex=legendCex, box.lwd=0, bty='n')

     legend("bottomright"
            , legend=c(agwGlue("(SL > 0) baseline: ",numSensitiveOverall,"/",numValidDiffPoints,": ",round(baselinePercentRight,0),"%")
              , agwGlue("Diamonds = In module. Recovery: ",numSensitiveInModule,"/",numValidInModule,": ",round(inModulePercentRight,0),"%")
              , agwGlue("Triangle = Top ",numCheckedVec[1],". Recovery: ",numCorrectVec[1],"/",numCheckedVec[1],": ",round(recoveredPercentVec[1],0),"%"
                        , " (",round(recoveredPercentVec[1]-baselinePercentRight,0),"%)"
                        , ", (",round(recoveredPercentVec[1]-inModulePercentRight,0),"%)"
                        )
              , agwGlue("Triangle = Top ",numCheckedVec[2],". Recovery: ",numCorrectVec[2],"/",numCheckedVec[2],": ",round(recoveredPercentVec[2],0),"%"
                        , " (",round(recoveredPercentVec[2]-baselinePercentRight,0),"%)"
                        , ", (",round(recoveredPercentVec[2]-inModulePercentRight,0),"%)"
                        )
              )
            , fill=c("white"
              , inModulePtFillColor
              , highScoringColorsVec[1]
              , highScoringColorsVec[2])
            , cex=legendCex, box.lwd=0, bty='n')
     
     ## Plot the "is differentially sensitive" results ===========
     if (yesSuper) { superSmoothed <- supsmu(slOverVec, isSensVec) }
     if (yesLoess) { loeSmoothed <- loess.smooth(slOverVec, isSensVec, family="gaussian", evaluation=80) }
     agwSmoothedNarrow <- agwSlidingWindowSmoother(slOverVec, isSensVec, windowSize=agwSmoothWindowSize, numSteps=agwSmoothNumSteps)

     agwWideWindowSize <- agwMax(2*agwSmoothWindowSize, ((agwMax(slOverVec)-agwMin(slOverVec))/10.0)) ## window size = 1/10 of the total range
     agwSmoothedWide   <- agwSlidingWindowSmoother(slOverVec, isSensVec, windowSize=agwWideWindowSize, numSteps=agwSmoothNumSteps)

     agwRawUnsmoothedData <- list(x=slOverVec, y=isSensVec, n=rep(1, times=length(slOverVec))) # typically not all that useful, but here's the raw data saved for every pathway!

     agwSmoothedAccum  <- agwSlidingWindowSmootherCumulative(slOverVec, isSensVec) # <-- note: does not make use of a window size at all!

     agwPrint("Loading module with: ", moduleName)

     gvCompleteCollection[[moduleName]] <<- list(           moduleName = moduleName
                                                 ,     slBaselineFract = baselinePercentRight/100
                                                 , moduleBaselineFract = inModulePercentRight/100
                                                 ,     cumulativeRecov = agwSmoothedAccum
                                                 ,      smoothedNarrow = agwSmoothedNarrow
                                                 ,        smoothedWide = agwSmoothedWide
                                                 ,       rawUnsmoothed = agwRawUnsmoothedData
                                                 ,            numSteps = agwSmoothNumSteps
                                                 ,        narrowWindow = agwSmoothWindowSize
                                                 ,          wideWindow = agwWideWindowSize)
     
     ## Complete collection: save the relevant statistics about this plot
     ## We're going to over-plot multiple charts later
     agwPrint("gvCompleteCollection length is now: ", length(gvCompleteCollection), " after adding module <", moduleName, ">.")
     
     screen(2) ; plot(c(-100), c(-100), type='n'
                      , xlim=plotXlim
                      , ylim=c(-1,1)
                      , xlab=xHypergeoLabel
                      , ylab="Fraction of diff. sens. ORF mutants"
                      , main=agwGlue("Fraction of ORFs that were diff. sensitive\n(i.e., (med(log10(growth_in_hitting_drugs)) - med(log10(growth_elsewhere))))/(med.abs.deviation) > ",round(sensThresh,1),")"))
     
     rect(xleft=0, ybottom=-1, xright=uninterestingCutoff+0.1, ytop=1, col=uninterestingOverlapColor, border=0) ## Fill the "uninteresting" areas
     rect(xleft=uninterestingCutoff, ybottom=-1, xright=interestingCutoff, ytop=1, col=slightlyInterestingOverlapColor, border=0)

     abline(h=0,col="black",lwd=zeroLineLwd)
     abline(h=c(0.25, 0.50, 0.75),col="gray",lwd=c(thinDelimiterLineLwd,thinDelimiterLineLwd,thinDelimiterLineLwd), lty=c(thinDelimiterLineStyle, thinDelimiterLineStyle, thinDelimiterLineStyle))

     points(slOverVec, slDiffVec/(agwMax(plotYlim,-plotYlim)), col=isSensColorVec, pch=isSensShapeVec, cex=pointScalingFactor, bg=isSensColorVec)


     if (yesLoess) { lines(loeSmoothed,   col=loeSmoothedColor,    type='o', cex=linePtCex, pch=PCH.SMALL.CIRCLE, lwd=smoothedLineWidth, lty="solid") } ## smoothed line
     if (yesSuper) { lines(superSmoothed, col=superSmoothedColor,  type='o', cex=linePtCex, pch=PCH.SMALL.CIRCLE, lwd=smoothedLineWidth, lty="solid") } ## smoothed line
     ##lines(agwSmoothedNarrow$x, agwSmoothedNarrow$y,           type='o', cex=linePtCex, pch=PCH.SMALL.CIRCLE,col=agwSmoothNarrowColor, lwd=smoothedLineWidth, lty="solid")
     lines(agwSmoothedWide$x, agwSmoothedWide$y,   type='o', cex=linePtCex,pch=PCH.SMALL.CIRCLE,col=agwSmoothWideColor, lwd=smoothedLineWidth)

     ##lines(agwSmoothedAccum$x, (agwSmoothedAccum$y - (baselinePercentRight/100)), type='o', cex=linePtCex,pch=PCH.SMALL.CIRCLE,col=agwAccumBaseColor, lwd=smoothedLineWidth)
     
     lines(agwSmoothedAccum$x, agwSmoothedAccum$y, type='o', cex=linePtCex,pch=PCH.SMALL.CIRCLE,col=agwSmoothAccumColor, lwd=smoothedLineWidth)

     ## Show the worse-than-baseline region--how well we do if we don't use SL data at all
     rect(xleft=0, ybottom=-100, xright=100, ytop=(inModulePercentRight/100), col=baselineModuleColor, border=0)  ## Plot the "baseline" shaded region
     ##abline(h=(inModulePercentRight/100),col=baselineModuleColor,lwd=baselineLwd, lty=baselineLty)

     ##abline(h=(baselinePercentRight/100),col=baselineSLColor,lwd=baselineLwd, lty=baselineLty)

     bottomRightLegendTextVec = c(
     agwGlue("In-module baseline (",round(inModulePercentRight,0),"%)")
     ##, agwGlue("Box-smoothed, width=",agwSmoothWindowSize)
     , agwGlue("Box-smoothed, width=", round(agwWideWindowSize,1))
     , agwGlue("Cumulative recovery rate"))
     ##, agwGlue("Cumulative - module baseline"))
     if (yesSuper) { bottomRightLegendTextVec = c(bottomRightLegendTextVec, "\"Super\" smoothed line") }
     if (yesLoess) { bottomRightLegendTextVec = c(bottomRightLegendTextVec, "Loess-smoothed line") }

     bottomRightLegendColorVec = c(baselineModuleColor, agwSmoothWideColor, agwSmoothAccumColor) ##, agwAccumBaseColor)
     if (yesSuper) { bottomRightLegendColorVec = c(bottomRightLegendColorVec, superSmoothedColor) }
     if (yesLoess) { bottomRightLegendColorVec = c(bottomRightLegendColorVec, loeSmoothedColor) }

     legend("bottomright"
            , legend=bottomRightLegendTextVec
            , fill=bottomRightLegendColorVec
            , cex=legendCex, title=NULL, box.lwd=0, bty='n')
     #points(slOverVec, isSensVec, col="blue",pch='X')
     
     ## Done plotting the "is differentially sensitive" results
     
     close.screen(all=TRUE)
     agwFinishPlot()

     ## Violin Plot: require(UsingR)
     ## simple.violinplot(decrease ~ treatment, data = OrchardSprays,col="bisque",border="black")

     ##  tm <- "OUR_KNOWN_TARGET:HAND_CURATED:microtubules"
     ##  gvModVec[[tm]]

     ##  colnames(gvGeneSLPathway)
     ##[1] "Deletion_Gene"    "Module"           "SL_Overlap"
     ##[4] "N_Intersect"      "N_Gene_SL"        "N_Module_Members"
     ##[7] "N_Universe"
     
     ## This is the X axis, indicating the overlap with modules
     ##overlaps <- gvGeneSLPathway[gvModVec[[tm]][1:length(gvModVec[[tm]])],"SL_Overlap"]
     ##gnam <- gvGeneSLPathway[gvModVec[[tm]][1:length(gvModVec[[tm]])],"Deletion_Gene"]
     ##plot(density(overlaps))
     ##plot(x=overlaps, y=overlaps, main="About the <MODULE>: Each point is an ORF", xlab="Hypergeometric score for the ORF's SGI overlap with the module in question")
     ##text(x=overlaps[1:10], y=overlaps[1:10], gnam[1:10],pos=POS.RIGHT)
}


agwFigureA_DetailPlots_main <- function() { ## <-- the name here is important!
     agwFigureA_GenerateManyDetailPlots() ## <-- call the function!!
}




##
