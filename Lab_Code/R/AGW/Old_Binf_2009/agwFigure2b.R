
## To run:         source("/Users/alexgw/R/src/agwFigure2b.R")      ;

## This makes the big histogram of drug differential sensitivities across modules.

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

FINALIZE = TRUE

shortNamesForModulesHash <- agwNewHash()
agwHashPutMultiple(shortNamesForModulesHash
                   , keys = c("OUR_KNOWN_TARGET:HAND_CURATED:AGW_Target_from_Imming_DNA_RNA_synth"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:amsacrine"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:microtubules"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:bleomycin"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:camptothecin"
                     , "GoDDCR"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:cell_wall"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:lipid_biosynthesis"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:AGW_Target_from_Imming_sterol_activity"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:peptidyl-prolyl_cis-trans isomerase activity"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:floxuridine"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:actin"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:AGW_Target_from_Imming_oxidoreductase_regulation"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:methoxsalen"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:sphingolipid_biosynthesis"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:TOR_signaling"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:AGW_Target_from_Imming_eletron_transport"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:chromatin"
                     , "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:staurosporine"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:AGW_Target_from_Imming_actin"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:kinase_inhibitor"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:oxidation"
                     , "OUR_KNOWN_TARGET:HAND_CURATED:phosphatase_inhibitor"
                     )
                   , values = c("DNA / RNA Synthesis (from Imming)"
                     , "Target of Amsacrine (Drugbank)"
                     , "Microtubules"
                     , "Target of Bleomycin (Drugbank)"
                     , "Target of Camptothecin (Drugbank)"
                     , "DNA Damage Checkpoint Response (DDCR)"
                     , "Cell Wall"
                     , "Lipid Biosynthesis"
                     , "Sterol Activity (from Imming)"
                     , "Peptidyl-prolyl Cis-trans Isomerase Activity"
                     , "Target of Floxuridine (Drugbank)"
                     , "Actin"
                     , "Oxidoreductase Regulation (from Imming)"
                     , "Target of Methoxsalen (Drugbank)"
                     , "Sphingolipid Biosynthesis"
                     , "TOR Signaling"
                     , "Electron Transport"
                     , "Chromatin"
                     , "Target of Staurosporine (Drugbank)"
                     , "Actin (from Imming)"
                     , "Kinase Inhibitor"
                     , "Oxidation"
                     , "Phosphatase Inhibitor"
                     )
                   )


agwFigure2bBoxplot <- function(plotType="boxonly", drawBackground=TRUE) {

     stopifnot(plotType == "violin" || plotType == "boxonly"
               || plotType == "boxpoints" || plotType == "pointsonly")
     
     agwGlobalLoad("gvNeighborStatsList")
     agwGlobalLoad("gvPathwayMembershipHash")
     boxWidthWex      <- NA
     stripchartOffset <- NA
     stripchartCex    <- NA
     stripchartJitter <- NA
     boxPch           <- NA
     missedThreshColorTransparent <-  "#00000099"
     madeThreshColorTransparent   <-  "#FF000099"
     if (plotType == "boxonly") {
          boxWidthWex <- 0.7
          boxPch      <- PCH.SMALL.CIRCLE
     } else if (plotType == "violin") {
          boxPch = "Y"  #rnorm(5)
          boxWidthWex = 1
     } else if (plotType == "boxpoints") {
          boxWidthWex      <-  0.30
          stripchartOffset <- -0.25
          stripchartCex    <-  0.40
          stripchartJitter <-  0.03
          boxPch           <- PCH.SMALL.CIRCLE
     } else if (plotType == "pointsonly") {
          stripchartOffset   <- 0
          boxWidthWex        <- 0
          stripchartCex      <- 1.00
          stripchartJitter   <- 0.15
          boxPch             <- ' ' ## no plotting of the points in the boxplot part--only plot the points in "stripchart" separately
          missedThreshColorTransparent <- "#000000AA"
          madeThreshColorTransparent   <- "#FF0000FF"
     }
     
     
     calcDifferenceVec <- function() {
          ## Calculates a vector that tells us how to sort the barplots
          ## in order to get a good most-to-least-interesting-performance-by-PGI-score ordering.
          retValVec <- vector(length=length(gvNeighborStatsList)) ## difference in sensitivity between made-the-threshold and did-not-make-the-threshold items. Only used for sorting!
          for (i in 1:length(gvNeighborStatsList)) {
               targetingGeneList  <- gvNeighborStatsList[[i]]
               madeVec <-  as.vector(sapply(targetingGeneList, function(a) { a$overlap >= kPGI_CUTOFF }))  ## this is all the data for genes that made the PGI overlap threshold
               failVec <-  as.vector(sapply(targetingGeneList, function(a) { a$overlap < kPGI_CUTOFF  }))  ## this is all the data for genes that did not make the PGI overlap threshold
               failVsPassDifference            <- NA
               requiredNumDataPointsForSorting <- 2
               if ((sum(madeVec, na.rm=TRUE) >= requiredNumDataPointsForSorting)   &&  (length(na.omit(targetingGeneList[madeVec])) >= requiredNumDataPointsForSorting)) {
                    cutoffPassedVec <- calcManyDiffSens(targetingGeneList[madeVec])
                    cutoffFailedVec <- calcManyDiffSens(targetingGeneList[failVec])
                    ##failVsPassDifference <- t.test(sapply(targetingGeneList[madeVec], function(a) { , targetingGeneList[failVec])
                    failVsPassDifference  <- mean(cutoffPassedVec, na.rm=TRUE) - mean(cutoffFailedVec, na.rm=TRUE)
                    
               } else {
                    ## not enough data points
               }
               retValVec[i] <- failVsPassDifference ## Used later for sorting!
          }
          return(retValVec) ## Used to sort the barplots in order of most-to-least difference between means
     }
     
     kPGI_CUTOFF <- 3.0 ## -log10(p-value) --- so 3 means p=0.001
     
     numBarsInGroup    <- 8 ## <-- set this manually!
     zeroThreshColor   <- "black" ;
     missedThreshColor <- "blue" ;
     madeThreshColor   <- "red"   ;
     
     inModuleColor    <- "#FFFF00" # yellow
     outModuleColor   <- "#4444FF" # blue
     anyModuleColor   <- "#22FF22" # green
     zeroScoreBGColor <- "#99999999" ;
     ishScoreBGColor  <- "#55559977" # these are scores between 0 and the PGI cutoff. They are nonzero-ish but not great. Hence "ish"
     goodScoreBGColor <- "#99FF9999" 
     barColors           <- c(outModuleColor,   inModuleColor
                              , outModuleColor, inModuleColor
                              , outModuleColor, inModuleColor
                              , "green"       , "green")     #,   anyModuleColor,   anyModuleColor)
     barBorderColors     <- c(zeroThreshColor, zeroThreshColor
                              , missedThreshColor, missedThreshColor
                              , madeThreshColor, madeThreshColor
                              , zeroThreshColor, madeThreshColor) # the color of the bars
     barBackgroundColors <- c(zeroScoreBGColor, zeroScoreBGColor
                              , ishScoreBGColor, ishScoreBGColor
                              , goodScoreBGColor, goodScoreBGColor
                              , "#FFFF3366", "#FF333366")

     stripchartPointColors <- barBorderColors

     stopifnot(length(barColors)       == numBarsInGroup)
     stopifnot(length(barBorderColors) == numBarsInGroup)
     
     theDiffVecHolder <- list()
     
     numUniqueDrugsWithDataToPlot  <- 0
     namesToUseVec    <- vector()
     boxplotOrderVec  <-  calcDifferenceVec() ## difference in sensitivity between made-the-threshold and did-not-make-the-threshold items. Only used for sorting the boxplots vertically! Not used within any *individual* plot
     orderedStatsList <- gvNeighborStatsList[ order(boxplotOrderVec) ] ## reorder the pathways based on the highest difference between high-PGI and low-PGI differential sensitivity
     rm(boxplotOrderVec)
     
     for (i in 1:length(orderedStatsList)) {
          moduleName       <- names(orderedStatsList)[i]
          longAbridgedName <- agwAbridgeModuleNames(moduleName)
          veryShortName    <- agwHashGet(shortNamesForModulesHash, key=moduleName, notFoundValue=moduleName)
          pathwayOverlappingGeneList  <- orderedStatsList[[i]]
          madeVec <- as.vector(sapply(pathwayOverlappingGeneList, function(a) { a$overlap >= kPGI_CUTOFF }))  ## this is all the data for genes that made the PGI overlap threshold
          lowVec  <- as.vector(sapply(pathwayOverlappingGeneList, function(a) { (a$overlap < kPGI_CUTOFF) && (a$overlap > 0.0) }))  ## this is all the data for genes that were under the PGI overlap threshold, but didn't have an PGI score of 0 
          zeroVec <- as.vector(sapply(pathwayOverlappingGeneList, function(a) { is.na(a$overlap) || a$overlap == 0 }))  ## this is all the data for genes that had no PGI connectivity
          
          ## Note: the variable elsewhere in this file called "numBarsInGroup" must be set to whatever the number of lines is here (as in, the number of drugs/variants we do for each drug). So if you have an entry in theDiffVecHolder[["someDrug"]] and theDiffVecHolder[["someDrug-that-made-cutoff"]], then numBarsInGroup would be set to 2
          orfsInThisModuleList <- agwHashGet(hash=gvPathwayMembershipHash, key=moduleName)
          orfsInModuleHash     <- agwNewHash() ; agwHashPutMembership(orfsInModuleHash, keys=orfsInThisModuleList, setValue=TRUE)
          rm(orfsInThisModuleList)
          inModuleVec <- as.vector(sapply(pathwayOverlappingGeneList, function(a) { agwHashContains(orfsInModuleHash, a$orf) } ) )
          
          stopifnot(length(madeVec) == length(pathwayOverlappingGeneList)) # experiments with PGI score >= kPGI_CUTOFF
          stopifnot(length(lowVec)  == length(pathwayOverlappingGeneList))  # experiments with PGI score > 0 but < kPGI_CUTOFF
          stopifnot(length(zeroVec) == length(pathwayOverlappingGeneList)) # experiments with PGI score == 0
          stopifnot(sum(madeVec, lowVec, zeroVec) == length(pathwayOverlappingGeneList)) # each score should be true SOMEwhere in here
          stopifnot(length(inModuleVec) == length(pathwayOverlappingGeneList))
          
          requiredNumDataPointsMadeThreshold <- 2  ## At least this many PGI points
          
          if ((sum(madeVec, na.rm=TRUE) >= requiredNumDataPointsMadeThreshold)
              &&  (length(na.omit(pathwayOverlappingGeneList[madeVec])) >= requiredNumDataPointsMadeThreshold)) {
               madeCutoffValuesVec <- calcManyDiffSens(pathwayOverlappingGeneList[madeVec])
               lowCutoffValuesVec  <- calcManyDiffSens(pathwayOverlappingGeneList[lowVec])
               zeroCutoffValuesVec <- calcManyDiffSens(pathwayOverlappingGeneList[zeroVec])
               
               ## In R, lists are populated in the order that their items are defined.
               ## So each one of these appends an item to the end of a list.
               ## Keep in mind that this unfortunately means that when plotting, the BOTTOM item shows up at the TOP of the plot.
               ## So these appear in reverse order (top to bottom -> bottom to top) in the final figure
               theDiffVecHolder[[agwGlue(moduleName,"_OUT_ZERO")]] <- calcManyDiffSens(pathwayOverlappingGeneList[zeroVec & (!inModuleVec)])
               theDiffVecHolder[[agwGlue(moduleName,"_IN_ZERO")]]  <- calcManyDiffSens(pathwayOverlappingGeneList[zeroVec & inModuleVec])
               theDiffVecHolder[[agwGlue(moduleName,"_OUT_LOW")]]  <- calcManyDiffSens(pathwayOverlappingGeneList[lowVec & (!inModuleVec)])
               theDiffVecHolder[[agwGlue(moduleName,"_IN_LOW")]]   <- calcManyDiffSens(pathwayOverlappingGeneList[lowVec & inModuleVec])
               theDiffVecHolder[[agwGlue(moduleName,"_OUT_MADE")]] <- calcManyDiffSens(pathwayOverlappingGeneList[madeVec & (!inModuleVec)])
               theDiffVecHolder[[agwGlue(moduleName,"_IN_MADE")]]  <- calcManyDiffSens(pathwayOverlappingGeneList[madeVec & inModuleVec])
               theDiffVecHolder[[agwGlue(moduleName,"_ANY_NOT_MADE")]] <- calcManyDiffSens(pathwayOverlappingGeneList[!madeVec])
               theDiffVecHolder[[agwGlue(moduleName,"_ANY_MADE")]]     <- calcManyDiffSens(pathwayOverlappingGeneList[madeVec])
               
               x <- theDiffVecHolder[[agwGlue(moduleName,"_ANY_MADE")]]
               y <- theDiffVecHolder[[agwGlue(moduleName,"_ANY_NOT_MADE")]]
               tCalc <- list("p.value"=1.0) ## worst possible score...
               #browser()
               if (length(x) > 2 && length(y) > 2) {
                    tCalc <- t.test(x, y, alternative="greater", na.action=na.exclude)
                    agwPrint("T-test p-value for differential sensitivities for PGI >= cutoff vs. diff. sens. for PGI < cutoff: ", tCalc$p.value)
               }
               
               print(tCalc$p.value)
               
               ##      theDiffVecHolder[[agwGlue(moduleName,"_ANY_FAILED")]] <- failedCutoffValuesVec ## failed the cutoff!
               ##      theDiffVecHolder[[agwGlue(moduleName,"_ANY_MADE")]] <- madeCutoffValuesVec   ## made the cutoff!
               ##      theDiffVecHolder[[agwGlue(moduleName," all MODULE members ")]] <- calcManyDiffSens(pathwayOverlappingGeneList[inModuleVec]) # failed the cutoff, AND in the module
               baseIndex   <- (numBarsInGroup*(numUniqueDrugsWithDataToPlot)) + 1
               scoreString <- format(kPGI_CUTOFF, nsmall=1)
               ## Note: these have to come in the same order as the definitions in "theDiffVecHolder" above
               beforeSize <- length(namesToUseVec)
               
               moduleAndDetails <- agwGlue(veryShortName, " (", length(orfsInModuleHash), " ORFs)")
               namesToUseVec[baseIndex+7] <- agwGlue(moduleAndDetails, "\n"
                                                     , "T-test p-value for PGI connectivity increasing sensitivity: ", round(tCalc$p.value, digits=4)
                                                     , "\nAny ORF with PGI score >= ", scoreString
                                                     , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_ANY_MADE")]])), ")")
               
               #print(paste("About to use..." , format(tCalc$p.value, nsmall=4)))
               namesToUseVec[baseIndex+6] <- agwGlue("Any ORF with PGI score < ", scoreString
                                                     , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_ANY_NOT_MADE")]])), ")")

               if (!FINALIZE) {
                    namesToUseVec[baseIndex+5] <- agwGlue(moduleAndDetails, "\n"
                                                          , "In-module, score >= ", scoreString
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_IN_MADE")]])), ")")
                    
                    namesToUseVec[baseIndex+4] <- agwGlue("Non-module, score >= ", scoreString
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_OUT_MADE")]])), ")")
                    
                    namesToUseVec[baseIndex+3] <- agwGlue(moduleAndDetails, "\n"
                                                          , "In-module, score between 0 and ",  scoreString
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_IN_LOW")]])), ")")
                    
                    namesToUseVec[baseIndex+2] <- agwGlue("Non-module, score between 0 and ", scoreString
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_OUT_LOW")]])), ")")
                    
                    namesToUseVec[baseIndex+1] <- agwGlue(moduleAndDetails, "\n"
                                                          , "In-module, zero PGI score "
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_IN_ZERO")]])), ")")
                    
                    namesToUseVec[baseIndex+0] <- agwGlue("Non-module, zero PGI score "
                                                          , " (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_OUT_ZERO")]])), ")")
               }
               ##      namesToUseVec[baseIndex+4] <- agwGlue(veryShortName, ", < ", scoreString, "  (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_2")]])), ")")
               ##      namesToUseVec[baseIndex+5] <- agwGlue(veryShortName, ", >= ", scoreString, "  (N=", length(na.omit(theDiffVecHolder[[agwGlue(moduleName,"_1")]])), ")")
               ##namesToUseVec[baseIndex+4] <- agwGlue("In module, any score")
               stopifnot((beforeSize + numBarsInGroup) == length(namesToUseVec)) ## numBarsInGroup must be the number of bars we just added here!
               
               numUniqueDrugsWithDataToPlot <- (1+numUniqueDrugsWithDataToPlot)
               agwPrint("Adding ", moduleName, " (fail-vs-pass PGI cutoff difference = ", kPGI_CUTOFF, ") to the boxplot...")
          } else {
               agwPrint("Skipping ", moduleName, " as it did not have any ORF mutants with an adequate PGI score (Not enough data points).")
          }
     }
     ### ============================
     
     stopifnot(length(theDiffVecHolder) == length(namesToUseVec))
     
     theRes <- 110           ## Resolution for the plot. Affects relative text size
     numSpacersBetweenGroups = 1  ## Number of spacers between bar graphs
     
     boxLocations <- vector()
     ##print("DIFF:") ; print(length(theDiffVecHolder)) ; print("DRUGS") ; print(numUniqueDrugsWithDataToPlot)
     stopifnot((numBarsInGroup*numUniqueDrugsWithDataToPlot) == length(theDiffVecHolder)) # <-- see above for notes on setting numBarsInGroup
     for (i in 1:numUniqueDrugsWithDataToPlot) {
          startLoc = (i-1)*(numBarsInGroup+numSpacersBetweenGroups)+1
          endLoc   = startLoc + numBarsInGroup - 1
          boxLocations <- c(boxLocations, seq(startLoc, endLoc, by=1))
     }
     
     boxLocations <- (1 - (numSpacersBetweenGroups/(numBarsInGroup+numSpacersBetweenGroups))) * boxLocations ## Horizontal boxplots have a bug! They don't properly allow non-regular "put a box at..." spacing. So we need to crunch them down a bit.
     #boxLocations <- seq(1:length(boxLocations)) #length(boxLocations)) ;
     
     bgString <- "";
     if (!drawBackground) {
          bgString <- "_no_bkground"
     }
     
     agwPreparePlot(directory=agwGlobalOutPath("Figure2b")
                    , file=agwGlue("Figure2b_Boxplot_", plotType, bgString)
                    , res=theRes, width=3000, height=160*length(theDiffVecHolder))
     m.top = (theRes/10); m.left = (theRes/3); m.bottom = (theRes/10); m.right = (theRes/20)
     par(mar=c(m.bottom, m.left, m.top, m.right))
     
     boxy <- function(add=FALSE) { ## <-- run the boxplot!
          boxplot(theDiffVecHolder
                  , names=namesToUseVec
                  , horizontal=TRUE
                  , main="Figure 2b: Differential Sensitivity by module\n(Sorted by the difference in differential sensitivity\nbetween high-scoring and low-scoring ORFs)"
                  , lwd=1
                  , boxwex=boxWidthWex ## box width amount (not a factor!)
                  , col=barColors
                  , pch=boxPch
                  , cex=0.8
                  , las=2 ## <-- perpendicular axis labels
                  , border=barBorderColors
                  , ylab=""
                  , xlab="Differential Sensitivity (difference in medians). Higher = \"actually targeting\" drugs are\nmore lethal. Growth is on a log_10 scale.)"
                  , at=boxLocations
                  , varwidth=FALSE
                  , lwd=4
                  , staplewex=1.0
                  , add=add)
     }
     
     if (plotType == "violin") {
          require(lattice)
          vioData <- list("data" = vector(), "group" = vector())
          
          kNUM_PADDING = ceiling(log10(length(theDiffVecHolder) + 5)) # how many digits?
          for (a in 1:length(theDiffVecHolder)) {
               blen <- length(theDiffVecHolder[[a]])
               theName <- agwGlue(format(a, width=kNUM_PADDING), " -- ", names(theDiffVecHolder)[a]) #originalModulePlotNames[a]
               ## above: Pad out the name for sorting. Sorting is done alphabetically, so numbers won't work (as "9" comes after both "10" and "1")
               theData <- theDiffVecHolder[[a]][1:blen]
               vioData$data  <- c(vioData$data, theData)
               vioData$group <- c(vioData$group, rep(theName, blen))
          }
          z <- bwplot(group ~ data, vioData  ## requires "lattice" package
                      , panel = function(..., box.ratio) {
                           panel.stripplot(..., jitter.data=TRUE, pch=PCH.X, cex=0.6, col="red")
                           panel.violin(...,   col = "#FFFF0044" #"transparent",
                                        , varwidth = FALSE, box.ratio = box.ratio)
                           ##panel.bwplot(...,  fill = "red", box.ratio = .1)
                           ##panel.fill(..., col=theV[index %% 2 + 1])
                      })
          ## plot.new()
          print(z)## <- this results in PLOTTING, but is NOT a standard plot call!! ## <-- have to put this here, or lattice won't plot see http://www.nabble.com/XYplot-in-Lattice-Package-td21491296.html
          ## "Lattice functions such as xyplot() create a graph object, but do not display it"
     }
     
     if (plotType != "violin") {
          boxy()
          barHeight <- 0.8
          if (drawBackground) {
               rect(xleft=rep(-100, length(theDiffVecHolder)), xright=rep(100,length(theDiffVecHolder))
                    , ybottom=boxLocations-(barHeight/2)
                    ,    ytop=boxLocations+(barHeight/2)
                    , col=barBackgroundColors
                    , border=NA)
          }
          abline(v=c(-2,2), lwd=1, col="#00000088", lty="dotted")
          #abline(v=0, lwd=2, col="orange", lty="solid")
          abline(v=0, lwd=2, col="black", lty="dashed")
          boxy(TRUE) ## plot the boxes a second time, on top of the background colors
     }
     
     
     if (plotType == "boxpoints" || plotType == "pointsonly") {
          stripchart(theDiffVecHolder, at=boxLocations+stripchartOffset
                     , method="jitter", jitter=stripchartJitter, vertical=FALSE, add=TRUE, pch=PCH.SMALL.CIRCLE, cex=stripchartCex
                     , col=stripchartPointColors)
     }
     
     legendCex <- 0.8
     
     if (plotType == "boxonly" || plotType == "boxpoints") {
          legend("topleft",
                 legend=c("Data from ORFs in each module"
                 , "Data from non-module ORFs")
                 , fill=c(inModuleColor, outModuleColor)
                 , cex=legendCex)
     }
     
     if (plotType == "boxonly" || plotType == "boxpoints" || plotType == "pointsonly") {
          legend("topright",
                 legend=c(agwGlue("ORFs with PGI score of at least ", scoreString)
                 , agwGlue("ORFs with PGI score less than ", scoreString))
                 , fill=c(goodScoreBGColor, ishScoreBGColor)
                 , cex=legendCex)
          
          legend("bottom",
                 legend=c("x-axis is in log_10 space.","x = 0 indicates equal median growth in both drug sets.","x = 0.3 indicates half median growth for mutants","exposed to the \"targeting\" drugs, relative","to median growth with other drugs.")
                 , cex=legendCex)
     }
     agwFinishPlot()
     ### ============================
}


agwFigure2b_main <- function() {
     ## Draw the images both with and without a background.
     
     for (b in c(TRUE, FALSE)) {
          agwFigure2bBoxplot("boxpoints", drawBackground=b)
          #agwFigure2bBoxplot("pointsonly", drawBackground=b)
          agwFigure2bBoxplot("boxonly", drawBackground=b)
     }
     agwFigure2bBoxplot("violin") ## Doesn't draw a background in any case, so the value doesn't matter.
}
