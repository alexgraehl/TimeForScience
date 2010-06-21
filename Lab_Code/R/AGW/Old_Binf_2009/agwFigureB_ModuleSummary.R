
##       source("/Users/alexgw/R/src/agwRankPlots.R");        ## <-- To reload this file

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");


RANK_PLOTS_DIRECTORY = "FigureB_ModSummary/"
kSGI_CUTOFF <- 3 ## -log10(p-value) --- so 3 means p=0.001.


#agwGlobalLoad("gvCompleteCollection") # <-- ??? maybe important


## =================================================
## Show the general informational plots in RANK_PLOTS_DIRECTORY

agwShowRankPlotsInVsOut <- function() {
     agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
     agwGlobalLoad("gvNeighborStatsList")
     agwGlobalLoad("gvNeighborStatsHash")
     #orfNameVec <- names(gvNeighborStatsList)
     moduleNameVec <- agwHashKeys(gvNeighborStatsHash)
     for (i in 1:length(moduleNameVec)) {
          #print("DEBUG: ONLY PRINTING ONE GRAPH!!!")
          moduleName      <- moduleNameVec[i]
          orfsWithDataVec <- agwHashKeys(agwHashGet(gvNeighborStatsHash, key=moduleName))
          statsList       <- sapply(orfsWithDataVec, function(a) { agwHashOfHashesGet(gvNeighborStatsHash, key1=moduleName, key2=a) })
          numGenes        <- length(orfsWithDataVec)
          hashExtractApply <- function(hash, theModuleName, inputKeys, listSubItem) {
               return (sapply(inputKeys, function(a) { agwHashOfHashesGet(hash=hash, key1=theModuleName, key2=a)[[listSubItem]]}));
          }
          sgiVec <- hashExtractApply(gvNeighborStatsHash, moduleName, orfsWithDataVec, "overlap") ## SGI overlap. Usually wilcoxon p-value
          
          inVec  <- hashExtractApply(gvNeighborStatsHash, moduleName, orfsWithDataVec, "rankMedIN")  ##  median RANK sensitivity of this orf among all the drugs that are said to target it. Lower means they actually were more-targeted
          outVec <- hashExtractApply(gvNeighborStatsHash, moduleName, orfsWithDataVec, "rankMedOUT")
          
          nDrugsTargetingModule       <- agwHashOfHashesGet(gvNeighborStatsHash, moduleName, orfsWithDataVec[1])$nTargetingModule
          nDrugsTargetingOtherModules <- agwHashOfHashesGet(gvNeighborStatsHash, moduleName, orfsWithDataVec[1])$nNotTargetingModule
          
          browser()
          
          ##t.test(inVec,outVec,paired=TRUE, alternative="greater") ## See if the PAIRED distributions really do favor the targeted-group being more sensitive than the non-targeted group
          
          nonNullInVec <- inVec[!sapply(inVec, is.null)]
          if ((numGenes > 0) && (length(nonNullInVec) > 0)) {
               worstPossibleRank <- max(as.vector(inVec,mode="double"), as.vector(outVec,mode="double"), na.rm=TRUE)
               agwPreparePlot(directory=agwGlue(agwGlobalOutPath(RANK_PLOTS_DIRECTORY)),file=agwGlue("Indicators_for_",moduleName,"_by_rank"), res=100, width=1500, height=1800)
               close.screen(all.screens=TRUE)
               split.screen(c(3,2))
               
               ## ==================================================
               screen(1)
               agwPrint("Generating a plot for ", moduleName, ". Red = orf sensitivity ranks for drugs that target the pathway, gray = orf sensitivity ranks among other drugs.")
               
               histLogMaxX  = 4
               histLogBreaks = seq(0,histLogMaxX, by=0.25)
               histLogXlim   = c(0, histLogMaxX)
               
               inColor  <- "#FF0000FF"
               outColor <- "#00005588"
               diffColor <- "#CC00FFFF"
               ##histLocal.xlim <- c(0, nrow(GLOBAL.HOM.SENS.COMPENDIUM))
               
               hist(log10(inVec)
                    , breaks=histLogBreaks, col=inColor, xlim=histLogXlim #histLocal.xlim
                    , xaxt="n"
                    , xlab="Sensitivity Rank (Rank 1 is most sensitive) (log scale)"
                    , border=0
                    , main=agwGlue("<", agwAbridgeModuleNames(moduleName), ">.\n Median rank is shown (1 = most sensitive to a drug).")
                    )
               
               hist(log10(outVec)
                    , breaks=histLogBreaks, col=outColor, xlim=histLogXlim #histLocal.xlim
                    , xaxt="n"
                    , border=0, add=TRUE)
               
               xMarksVec = seq(0, histLogMaxX, by = 1)
               axis(1, at=xMarksVec, labels=sapply(xMarksVec, function(a) { 10**a } ))
               legend("topleft"
                      , legend=c(
                        agwGlue("N = ", numGenes, " tested ORFs")
                        , "Median ORF sens. ranks in drugs TARGETING this pathway"
                        , "Median ORF sens. ranks in drugs NOT targeting this pathway"
                        ),
                      fill=c("black", inColor, outColor))
               ## ==================================================
               screen(2)
               sgiCutoff1 <- 2
               sgiCutoff2 <- 5
               
               uninterestingOverlapColor       <- "#DDEE99"
               slightlyInterestingOverlapColor <- "#EEFFCC"
               
               plotXLim <- c(0, agwMax(sgiVec))
               plotYLim <- c(worstPossibleRank,-worstPossibleRank)
               plot(sgiVec, inVec, pch=" ", xlim=plotXLim, ylim=plotYLim
                    , ylab="Rank (and Differential Rank)", xlab="SGI Centrality Score")
               
               rect(xleft=0, ybottom=agwMin(plotYLim), xright=sgiCutoff1, ytop=agwMax(plotYLim)
                    , col=uninterestingOverlapColor, border=0)
               rect(xleft=sgiCutoff1, ybottom=agwMin(plotYLim), xright=sgiCutoff2, ytop=agwMax(plotYLim)
                    , col=slightlyInterestingOverlapColor, border=0)
               
               markLWD = 5
               regularPointCex <- 1.0
               diffRankVec  <- (inVec - outVec) ## For each ORF, the median sensitivity *rank* to the drugs that target the pathway in question, minus the median sensitivity rank of the ORF to other drugs
               
               abline(h = 0, col="black", lwd=markLWD)
               
               superSmoothColor = "black"
               
               superRankSmoothed <- supsmu(sgiVec, diffRankVec)
               lines(superRankSmoothed, col=superSmoothColor,  type='o', cex=1.5, pch=PCH.SMALL.CIRCLE, lwd=3, lty="solid") ## smoothed line
               
               suppressWarnings(linearFit <- lm(diffRankVec ~ sgiVec))
               if (!is.null(linearFit)) {
                    abline(linearFit, col="black", lwd=2, lty="dashed")
               }
               
               points(sgiVec, inVec,  pch=PCH.DIAMOND, col=inColor, bg="gray", cex=regularPointCex)
               points(sgiVec, outVec, pch=PCH.DIAMOND,   col=outColor, bg="white", cex=regularPointCex)
               points(sgiVec, diffRankVec, pch=PCH.X, col=diffColor, bg="black",cex=1.1*regularPointCex)
               
               legend("topleft"
                      , legend=c("Med. rank among targeting drugs"
                      , "Med. rank among others"
                      , "(Med_IN - Med_OUT)"
                      , "Smoothed (Med_IN - OUT)")
                      , fill=c(inColor, outColor, diffColor, superSmoothColor)
                      )
               
               ## ========================================================================
               
               screen(3)
               sensList <- list(inVec, outVec)
               names(sensList) <- c(agwGlue("ORFs In targeting drugs (",nDrugsTargetingModule," drugs)")
                                    , agwGlue("ORFs In other drugs (",nDrugsTargetingOtherModules," drugs)")
                                    )
               boxplot(sensList, col=c(inColor, outColor)
                       , main=agwGlue(agwAbridgeModuleNames(moduleName))
                       , ylim=c(worstPossibleRank, 0)
                       , ylab="RANK (1 = most sensitive ORF)"
                       , xlab=agwGlue("N_ORFS = ", length(inVec))
                       )
               ## ========================================================================
               screen(4)
               plot(sgiVec,pch=PCH.SMALL.CIRCLE, ylab="Score", xlab="Index (in sorted vector)")
               
               ## ========================================================================
               screen(5)
               sortFromHighToLow <- TRUE
               rankedMatrix <- apply(GLOBAL.HOM.SENS.COMPENDIUM, BY.ROW, function(oneLineVec) {
                    cat("~") ;
                    z <- rank((-1 * oneLineVec), na.last=TRUE, ties.method="average"); ## <-- note, we NEGATE here, because we want to rank from high values (good ranks: 1,2,3,...) to low/negative values (bad ranks: 99,100,101,...). Also: NAs are NOT counted as equal for ranks!!!
                    numNonNA <- sum(!is.na(oneLineVec))
                    z[z > numNonNA] <- NA ## NA input data gets NA ranks
                    #print(numNonNA);
                    return(z)
               } )
               ##plot(sgiVec)
               agwFinishPlot()
               ## ========================================================================
          }
          ##par(ask=TRUE)
     }
}


## =======================================================================================
## =======================================================================================

agwFigureB_ModuleSummary_main <- function() {
     agwShowRankPlotsInVsOut()
}

#curve(sapply(gvNeighborStatsList[[1]],"[[", "rankMedOUT"), breaks=30, xlim=c(0,1000),col="red")


#testFF <- function() {
#  lapply(gvCompleteCollection, function(a) { agwPrint(a$cumulativeRecov$x, ", " , a$cumulativeRecov$y) })
#  
#}


##agwShowRankPlotsInVsOut()

