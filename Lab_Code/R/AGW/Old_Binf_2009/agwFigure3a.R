#
#
##       source("/Users/alexgw/R/src/agwFigure3a.R");        ## <-- To reload this file

## Here is an example where the "drugs that target the module" set (the "x" variable / first argument in t.test) causes significantly more lethality to the particular ORF deletion than "drugs that target other modules" (the "y" variable / second argument in t.test)
##ttt <- t.test(c(15, 15, 15, 15, 15), c(12,11,12,13,14,12,15,16,12,11,11,8,12,12,11,10),alternative="greater")
## Just for testing


#allORFs  <- rownames(gv.HOM.RANKS.BY.DRUG)
#allDrugs <- colnames(gv.HOM.RANKS.BY.DRUG)

##theIdx = 1

k_FIGURE_3_FOLDER_NAME <- "Figure3a"

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");


## =========================================================================


## =================================================

agwFigure3aPlot <- function() {
     kSGI_CUTOFF_FIG3a <- 3 ## -log10(p-value) --- so 3 means p=0.001
     agwGlobalLoad("gv.HOM.RANKS.BY.DRUG")
     agwGlobalLoad("gv.HOM.RANKS.BY.ORF")

     heatmapLegendCex <- 0.7
     heatmapLegendBG  <- "#FFFFFF99"

     agwGlobalLoad("gvNeighborStatsList")
     
     ## Figure 3A
     for (theIdx in 1:5) { #length(gvNeighborStatsList)) {
          moduleName    <- names(gvNeighborStatsList)[theIdx]
          print("DEBUGGGGGGGGGGGGGGGGGGGGGG");
          
          targetingDrugsVec <- names(gvDrugTargetMapping[[moduleName]]) ## drugs that are annotated as targeting this biochemical module
          agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
          testedDrugsVec          <- colnames(GLOBAL.HOM.SENS.COMPENDIUM)
          targetingTestedDrugsVec <- intersect(testedDrugsVec, targetingDrugsVec) ## drugs that target this module AND were tested
          
          theStatsList <- gvNeighborStatsList[[theIdx]]
          theStatsList <- theStatsList[1:min(length(theStatsList),250)] ### <- debugging!!!!
          print("DEBUGGGGGGGGGGGGGGGGGGGGGG:: LIMITING STATS LIST TO 250!!");
          
          orfsWithFullDataVec    <- intersect(colnames(gv.HOM.RANKS.BY.ORF), names(theStatsList))
          orfsInStatsListOnlyVec <- setdiff(names(theStatsList), colnames(gv.HOM.RANKS.BY.ORF))
          
          modDataList <- theStatsList[sort(orfsWithFullDataVec)] ## Note that this is ordered alphabetically!!
          
          if (length(orfsInStatsListOnlyVec) > 0) {
               print("We had to kick out any ORFs that we had data for in the gvNeighborStatsList but not in gv.HOM.RANKS.BY.ORF. Here is a list:")
               print(orfsInStatsListOnlyVec)
          }
          
          orderBySGIOverlapVec <- order(sapply(modDataList,"[[", "overlap")) ## sorted so that early values are NA and bad values, and later ones are good
          modDataList  <- modDataList[orderBySGIOverlapVec] ## now "modDataList" is sorted by SGI score
          sgiVec       <- sapply(modDataList,"[[", "overlap")
          geneNameVec  <- sapply(modDataList,"[[", "gene")
          orfNameVec   <- sapply(modDataList,"[[", "orf")
          
          agwGlobalLoad("gvPathwayMembershipHash")
          orfsInThisModuleVec         <- agwHashGet(hash=gvPathwayMembershipHash, key=moduleName)
          inModuleBoolVec <- (orfNameVec %in% orfsInThisModuleVec) # check to see which ORFs are in fact in the module
          
          sgiInModuleVec      <-      sgiVec[inModuleBoolVec] # <-- these are ONLY including the genes that are in the given GO category
          geneNameInModuleVec <- geneNameVec[inModuleBoolVec] # <-- these are ONLY including the genes that are in the given GO category
          orfNameInModuleVec  <-  orfNameVec[inModuleBoolVec] # <-- these are ONLY including the genes that are in the given GO category
          # The idea with the above "only in the GO category" vectors is to see what the performance of picking indicator genes is if you
          # *just* pick a gene in the module.
          
          ## =====================
          plotWidth <- 1000 + (length(geneNameVec)*12)
          agwPreparePlot(directory=agwGlobalOutPath(agwGlue(k_FIGURE_3_FOLDER_NAME))
                         ,    file=agwAbridgeModuleNames(moduleName)
                         , res=100, width=plotWidth, height=2500)
          
          screenSplitLocs <- c(0.85, 0.5, 0.35, 0.17, 0.10) ## the vertical values to split up the screen at
          close.screen(all=TRUE)
          split.screen( c(7, 1)) #matrix(data=c( ## left, right, bottom, top
                               #0.0, 1.0, screenSplitLocs[1], 1.00
                               #, 0.0, 1.0, screenSplitLocs[2], screenSplitLocs[1]
                               #, 0.0, 1.0, screenSplitLocs[3], screenSplitLocs[2]
                               #, 0.0, 1.0, screenSplitLocs[4], screenSplitLocs[3]
                               #, 0.0, 1.0, 0, screenSplitLocs[4]
                               #), byrow=T, ncol=4))
          
          stopifnot(all(orfNameVec == names(sgiVec))) # <-- make sure the names line up properly!
          
          ## ============================================
          m.top = 2; m.left = 8; m.bottom = 8; m.right = 2
          orfXLim <- c(-1, length(sgiVec)+1)
          orfx.cex.axis <- 0.8
          ## ============================================
          screen(1)
          xMarkVec <- seq(1,length(sgiVec))
          madeSGICutoffColor   <- "red"  # Colors on the top bar graph for the ORFs that made the SGI cutoff
          failedSGICutoffColor <- "dark green" # Colors on the top bar graph for the ORFs that did not make the SGI cutoff
          zeroSGICutoffColor <- "gray" # ORFs that failed the SGI cutoff, and in fact had ZERO overlap with the module
          sgiColorVec <- ifelse((sgiVec > kSGI_CUTOFF_FIG3a), madeSGICutoffColor, failedSGICutoffColor)
          sgiColorVec[sgiVec >= kSGI_CUTOFF_FIG3a] <- madeSGICutoffColor
          sgiColorVec[sgiVec < kSGI_CUTOFF_FIG3a]  <- failedSGICutoffColor
          sgiColorVec[sgiVec == 0]                 <- zeroSGICutoffColor
          par(mar=c(m.bottom, m.left, m.top, m.right))
          
          plotYLim = c(0, max(sgiVec)*1.05)
          plot(x=xMarkVec, y=sgiVec
               , xaxt="n"
               , xlab=""
               , ylab="SGI Overlap Score w/Module (log10(p))"
               , main=agwGlue("Figure 3A.  X-axis: ORFs with non-zero SGI overlap with ", agwAbridgeModuleNames(moduleName))
               , pch=PCH.SMALL.CIRCLE
               , col=sgiColorVec
               , xaxs="i"
               , xlim=orfXLim
               , ylim=plotYLim
               , type='h'
               , lwd=5
               )
          
          everyTenColor <- "#0000FF20"  ## the color for the iTunes-playlist-like boxes that visually separate every 10 ORFs
          boxStepSize <- 10
          b <- length(sgiVec)
          while (b > 0) { ## Draw iTunes-playlist-like boxes to visually separate every 10 ORFs
               rect(xleft=b-boxStepSize, xright=b, ybottom=0, ytop=plotYLim, col=everyTenColor, border=NA)
               b <- (b - (2*boxStepSize))
          }
          
          abline(h=kSGI_CUTOFF_FIG3a, lwd=2, col=madeSGICutoffColor, lty="dashed")
          par(las=2)
          axis(1, at=xMarkVec, labels=agwGlue("(", orfNameVec, ")  ", geneNameVec), cex.axis=orfx.cex.axis)
          
          ## ==================================
          screen(2)
          zMat <- gv.HOM.RANKS.BY.ORF[, which(colnames(gv.HOM.RANKS.BY.ORF) %in% orfNameVec)] # first: the drugs, second, the orfs. So it's [DRUG, ORF]
          
          allDrugsVec <- rownames(zMat)
          dTarVec     <- which(rownames(zMat) %in% targetingTestedDrugsVec)
          
          toSort <- matrix(data=rep(0,length(zMat)) ## <-- fill it with zeroes!
                           , nrow=nrow(zMat)
                           , ncol=ncol(zMat)
                           , dimnames=list(rownames(zMat),colnames(zMat))
                           )
          
          maxDrugRankPerORFVec <- apply(zMat, BY.COL, max, na.rm=TRUE) ## For each ORF, what's the very worst rank for any drug? This can be different for each ORF, due to ties and various things.
          
          accum <- 1 ## <-- start at 1! Zero means "no drug at all"
          drugColorsVec <- vector()
          for (d in dTarVec) {
               toSort[d,] <- sapply(toSort[d,],function(a){ return(accum) }) ## Then mark the drugs that ARE targeting this particular module
               accum <- (1+accum)
          }
          colorForTargetingDrugs <- "black"
          drugColorsVec <- rep(colorForTargetingDrugs, length(dTarVec)) #rainbow(n=length(dTarVec))
          
          tsMat <- apply(zMat, BY.COL, function(a) {
               theOrder <- order(a) ;
               return(toSort[theOrder]) } ) ## Sort each column so that the drugs are sorted from most-to-least sensitive (uses the *ranks* from "modDataList"--rank 1 = most sensitive)
          
          tsNumRows <- nrow(tsMat)
          for (ss in 1:ncol(tsMat)) { ## For each ORF...
               numToReplace <- abs(maxDrugRankPerORFVec[ss] - tsNumRows) + 1
               if (numToReplace > 0) {
                    tsMat[maxDrugRankPerORFVec[ss]:tsNumRows, ss] <- rep(-1, numToReplace); ## These are the "invalid" data regions that are worse than the maximum rank
               }
          }
          
          tsMat <- tsMat[, order(colnames(tsMat))] ## reorder columns alphabetically. We don't know how they were ordered before!
          tsMat <- tsMat[, orderBySGIOverlapVec] ## now the ORF columns are ordered by SGI score, so they match the same order seen in "sgiVec"
          
          stopifnot(NULL == rownames(tsMat)) ## tsMat does not have rownames! each column is sorted INDEPDENDENTLY, so it doesn't make any sense to say that there are rownames.
          worstRank <- max(zMat, na.rm=TRUE)
          stopifnot(0 == sum(is.na(tsMat))) ## there should be no NA values in tsMat
          finalMat <- t(tsMat)
          stopifnot(length(sgiVec) == nrow(finalMat))
          
          top50              <- 50    ##worstRank*(1.0 - 0.90)
          top50Color         <- "red"
          backgroundForCells <- "white" ## background for the "show where the various drugs fell" background
          par(mar=c(m.bottom, m.left, m.top, m.right))
          image(x=1:nrow(finalMat), y=1:ncol(finalMat), z=finalMat
                , main=agwGlue(length(targetingTestedDrugsVec), " tested drugs target the module \"", agwAbridgeModuleNames(moduleName), "\"")
                , xlab="", ylab=agwGlue("Drug Sens. Rank (lower = drug is more deadly). No. drugs: ", worstRank, ". Colored bars indicate drugs that target this module")
                , xaxt="n"
                , xlim=orfXLim
                , ylim=c(worstRank+1,1)
                , col=c(backgroundForCells)
                , lwd=5
                )
          ## browser()
          ## rect(xleft=orfXLim[1], xright=orfXLim[2], ytop=0, ybottom=top50, col=top50Color, border=NA)
          
          abline(h=top50, col=top50Color, lwd=3, lty="dashed") ## show a line marking off the top 50 drugs
          ##rect(xleft=orfXLim[1], xright=orfXLim[2], ytop=0, ybottom=percentile95, col="#FFFF99FF", border=NA)
          
          b <- length(sgiVec)
          while (b > 0) { ## Draw iTunes-playlist-like boxes to visually separate every 10 ORFs
               rect(xleft=b-boxStepSize, xright=b, ybottom=0, ytop=(worstRank+1), col=everyTenColor, border=NA) ;
               b <- (b - (2*boxStepSize))
          }
          
          noDataColor          <- "#00000000" # <-- this is (effectively) the background color for the entire chart. It can have transparency!
          beyondWorstRankColor <- "#AAAAAA" ## These are ranks that are IMPOSSIBLE for a drug to get, because of ties! (So if there are 1000 drugs, the worst rank could still be 500 if ties got the average rank and every single drug had the same value for an ORF)
          
          image(x=1:nrow(finalMat), y=1:ncol(finalMat), z=finalMat ## Plot the actual drugs
                , main=""
                , xlab="", ylab=""
                , xaxt="n"
                , yaxt="n"
                , col=c(beyondWorstRankColor, noDataColor, drugColorsVec)
                , cex.axis=orfx.cex.axis
                , add=TRUE
                )
          
          legend("topleft"
                 , legend=c("Top 50 most-sensitizing drugs")
                 ,   fill=c(top50Color)
                 ,     bg=heatmapLegendBG
                 ,    cex=heatmapLegendCex
                 )
          
          numDrugsTargetingThisModule <- length(targetingTestedDrugsVec)
          legend("bottomleft"
                 , legend=c(
                   agwGlue("Num ORFs in module: ", length(sgiVec))
                   , agwGlue("Num drugs targeting module: ", numDrugsTargetingThisModule)
                   , targetingTestedDrugsVec)
                 , fill=c("black", "black", drugColorsVec)
                 ,   bg=heatmapLegendBG
                 ,  cex=heatmapLegendCex
                 )
          ## axis(A.Y.AXIS, at=seq(1,ncol(zMat),by=y.label.every), labels=colnames(zMat)[seq(1,ncol(zMat),by=y.label.every)], las=1)
          
          orflabels <- rownames(finalMat)
          geneLabels <- sapply(orflabels, function(a) {
               return(agwGetAnnotationFromOrf(a)$gene)
          })
          
          axis(A.X.AXIS
               , at=seq(1,nrow(finalMat),by=1)
               , labels=agwGlue("(",orflabels, ") ", geneLabels)
               , cex.axis=orfx.cex.axis
               , las=2)
          
          ## =======================================================================================
          
          agwPlotNiftyMatrix <- function(theMat, mainLabel="", imageColorsVec, imageBreaksVec) {
               par(mar=c(m.bottom, m.left, m.top, m.right))
               image(x=1:nrow(theMat), y=1:ncol(theMat), z=theMat
                     , xlim=orfXLim
                     , ylim=c(ncol(theMat),1)
                     , main=mainLabel
                     , xlab=""
                     , ylab=""
                     , xaxt="n"
                     , yaxt="n"
                     , col=imageColorsVec
                     , breaks=imageBreaksVec
                )
               
               axis(A.Y.AXIS, at=seq(1,length(topKVec)),
                    , labels=agwGlue("Top ", topKVec)
                    , las=LAS.HORIZONTAL.TEXT)
               
               axis(A.X.AXIS, at=seq(1,nrow(finalMat))
                    , labels=agwGlue("(", orfNameVec, ")  ", geneNameVec)
                    , cex.axis=orfx.cex.axis
                    , las=2)
               
               legend("topleft"
                      , legend=agwImageBreakLabelsFromColors(imageColorsVec, imageBreaksVec)
                      , fill=imageColorsVec
                      , bg=heatmapLegendBG
                      , cex=heatmapLegendCex
                      )
               
          }
          
          screen(3)
          # Q: out of the top N drugs that are selected by this indicator, what is the total number of them that are actually supposedly targeting this pathway? (Raw numbers, not scores!)
          
          valuesGTEThanThisAreTargetingDrugs <- 1 ## -1 is an empty spot, and 0 is a drug that is *not* claimed to target this pathway. 1,2,3,... etc are the drugs that *do* target it.
          topKVec <- c(1,2,3,4,5,6,7,8,9,10,15,20,50)  ## How many top predicted genes do we want to look at? (i.e., run experiments and see if they are, in fact, good indicator genes)
          
          numInTopMat <- matrix(nrow=nrow(finalMat), ncol=length(topKVec)
                                , dimnames=list(rownames(finalMat), topKVec)) #paste("N in top", topKVec, "most sensitive drugs")))
          for (kIdx in 1:length(topKVec)) {
               k <- topKVec[kIdx]
               performVec <- apply(finalMat,BY.ROW,function(a) { sum( (a[1:k] > 0) ) }) ## figure out what number of the topK drugs are in fact targeting the desired pathway
               numInTopMat[,kIdx] <- performVec
          }
          
          numDivisions1    <- 9
          break1BreaksVec  <- c(-999, seq(0,numDivisions1-2), max(numInTopMat+1, numDivisions1-2)) ## -1 and 100 are guaranteed to be off the ends of the scale
          break1ColorsVec  <- agwGetHeatColors(length(break1BreaksVec)-1) #rev(gray((0:(length(break1BreaksVec)-2))/(length(break1BreaksVec)-2)))
          stopifnot(length(break1ColorsVec)+1 == length(break1BreaksVec))
          ## in numInTopMat: 0 means "didn't get any drugs in the topK"
          ## 1 means "got one drug", 2 is 2, etc...
          agwPlotNiftyMatrix(numInTopMat, mainLabel=agwGlue("Recovery of no. targeting drugs in top N (of ", length(targetingTestedDrugsVec), " possible drugs.")
                             , break1ColorsVec, break1BreaksVec)
          
          ## =======================================================================================
          screen(4)
          
          ## Out of the top N drugs selected by this indicator, what fraction are the known targets? (How precise are the top N?)
          indexInTopK <- 0
          precisionInTopMat <- apply(numInTopMat, BY.COL, function(a) { ## <-- a FRACTION, so it's between 0 and 1
               indexInTopK <<- (indexInTopK+1)
               return(a/topKVec[indexInTopK]) ## How many drugs we got, versus how many we checked in the top K (precision! out of K trials, we got N of them)
          })
          
          stopifnot(max(precisionInTopMat, na.rm=TRUE) <= 1) ; stopifnot(min(precisionInTopMat, na.rm=TRUE) >= 0)
          numDivisions2   <- 5
          #break2ColorsVec  <- heat.colors(numDivisions2+1)#c("black", "#772222", "red","orange","yellow","green","cyan","blue","#EE11BB","#FF66FF", "white")
          break2BreaksVec <- c(-999, 0, seq((1.0/numDivisions2),1.0,by=(1.0/numDivisions2))) ## -1 and 100 are guaranteed to be off the ends of the scale
          break2ColorsVec <- agwGetHeatColors(length(break2BreaksVec)-1) #rev(gray((0:(length(break2BreaksVec)-2))/(length(break2BreaksVec)-2)))
          stopifnot(length(break2ColorsVec)+1 == length(break2BreaksVec)) ## <-- required by "image"
          agwPlotNiftyMatrix(precisionInTopMat
                             , mainLabel=agwGlue("Precision: Out of the top N drugs, the fraction that were targeting-this-pathway drugs. (Out of ", length(targetingTestedDrugsVec), " drugs.)")
                             , break2ColorsVec
                             , break2BreaksVec)
          
          ## =======================================================================================
          screen(5)
          accumNumInTopMat <- numInTopMat
          for (topIdx in 1:ncol(numInTopMat)) {
               for (geneIdx in nrow(numInTopMat):1) {
                    numCheckedGenesByNow      <- nrow(numInTopMat) - geneIdx + 1  # how many gene deletion mutants have we checked by now?
                    numCheckedDrugsPerGene    <- topKVec[topIdx]
                    numberOfCheckedDrugsByNow <- numCheckedGenesByNow * numCheckedDrugsPerGene # how many top-most-sensitive drugs have we checked by now?
                    totalHittingDrugsInTop <- sum(numInTopMat[nrow(numInTopMat):geneIdx,topIdx]) # the number of hitting-the-pathway drugs in the top N-most-sensitive drugs. Note that "hitting" is an adjective here: "hitting drugs" are drugs that are annotated as hitting the module in question.
                    maxPossibleAcrossGenes <- numCheckedGenesByNow*min(numDrugsTargetingThisModule, topKVec[topIdx]) # how many top drugs COULD we have gotten
                    accumNumInTopMat[geneIdx, topIdx] <- totalHittingDrugsInTop
               }
          }
          
          accumDivisions <- numDrugsTargetingThisModule ## <-- this can be scaled up or down as well... nothing too special about this particular choice
          accumBreaksVec <- NA
          if (accumDivisions < 3) {
               accumBreaksVec <- c(-999, 0)
          } else {
               accumBreaksVec <- c(-999, seq(0, accumDivisions-3, by=1),max(accumNumInTopMat)+1) #seq(0, 1.0, by=(1.0/accumDivisions)) ## -1 and 100 are guaranteed to be off the ends of the scale
          }
          accumColorsVec <- agwGetHeatColors((length(accumBreaksVec)-1))
          
          agwPlotNiftyMatrix(accumNumInTopMat
                             , mainLabel=agwGlue("Cumulative recovery of no. targeting drugs in top N (of ", length(targetingTestedDrugsVec), " possible drugs)")
                             , imageColorsVec=accumColorsVec
                             , imageBreaksVec=accumBreaksVec)
          


          ## =======================================================================================
          screen(6) ## The "expected vs actual recover" graph
          
          expModuleMat <- matrix(nrow=length(sgiInModuleVec), ncol=length(topKVec)) # <-- "expected number of hitting drugs we find for N tested drugs"
          actModuleMat <- matrix(nrow=length(sgiInModuleVec), ncol=length(topKVec)) # <-- "actual number of hitting drugs we find for N tested drugs"
          maxModuleMat <- matrix(nrow=length(sgiInModuleVec), ncol=length(topKVec)) # <-- "max possible number of hitting drugs we find for N tested drugs"
          
          chanceOfGettingPositiveRandomly <- numDrugsTargetingThisModule/length(allDrugsVec) # chance of randomly selecting a drug and having it be a hitting drug
          

          #browser()
          calcVariousHittingDrugRecoveryMatrices <- function(sourceMat) {
               returnableExpMat <- matrix(nrow=nrow(sourceMat), ncol=ncol(sourceMat), dimnames=list(rownames(sourceMat), colnames(sourceMat))) # <-- EXPECTED number of hitting drugs we find for N tested drugs
               # ^ EXPECTED number of hitting drugs we find for N tested drugs
               returnableActMat <- matrix(nrow=nrow(sourceMat), ncol=ncol(sourceMat), dimnames=list(rownames(sourceMat), colnames(sourceMat))) # <-- ACTUAL number of hitting drugs found for testing N drugs
               # ^ ACTUAL number of hitting drugs found for testing N drugs
               returnableMaxMat <- matrix(nrow=nrow(sourceMat), ncol=ncol(sourceMat), dimnames=list(rownames(sourceMat), colnames(sourceMat))) # <-- MAX: the maximum (i.e., best) POSSIBLE number of recovered hitting drugs, in testing N drugs
               # ^ MAX: the maximum (i.e., best) POSSIBLE number of recovered hitting drugs, in testing N drugs

               returnableAccumMat <- matrix(nrow=nrow(sourceMat), ncol=ncol(sourceMat), dimnames=list(rownames(sourceMat), colnames(sourceMat)))
               for (i in 1:ncol(sourceMat)) { # i = top index
                    for (j in nrow(sourceMat):1) { # j = gene index
                         numCheckedGenesByNow      <- nrow(sourceMat) - j + 1  # how many gene deletion mutants have we checked by now?
                         numCheckedDrugsPerGene    <- topKVec[i]
                         numberOfCheckedDrugsByNow <- numCheckedGenesByNow * numCheckedDrugsPerGene # how many top-most-sensitive drugs have we checked by now?
                         totalHittingDrugsInTop    <- sum(sourceMat[nrow(sourceMat):j,i])
                         # Above: the number of hitting-the-pathway drugs in the top N-most-sensitive drugs. Note that "hitting" is an adjective here: "hitting drugs" are drugs that are annotated as hitting the module in question.
                         totalHittingDrugsInTopInModuleOnlyGenes <- sum(sourceMat[nrow(sourceMat):j,i])
                         maxPossibleAcrossGenes  <- numCheckedGenesByNow*min(numDrugsTargetingThisModule, topKVec[i]) # how many top drugs COULD we have gotten
                         numExpected             <- chanceOfGettingPositiveRandomly * numberOfCheckedDrugsByNow
                         returnableExpMat[j, i] <- numExpected            ## Number of drugs we would EXPECT to get...
                         returnableActMat[j, i] <- totalHittingDrugsInTop ## Actual number of drugs we got...
                         returnableMaxMat[j, i] <- maxPossibleAcrossGenes
                         returnableAccumMat[j, i] <- totalHittingDrugsInTop - numExpected #totalHittingDrugsInTop #/maxPossibleAcrossGenes #log(totalHittingDrugsInTop+0.1)
                    }
               }
               return(list(expected=returnableExpMat, actual=returnableActMat, max=returnableMaxMat, accum=returnableAccumMat))
          }
          
          recovery <- calcVariousHittingDrugRecoveryMatrices(numInTopMat)
          
          numInTopInModuleMat <- numInTopMat[orfNameInModuleVec,]
          recoveryInModule <- calcVariousHittingDrugRecoveryMatrices( numInTopInModuleMat )
          
          diffAccumBreaksVec <- c(-999, seq(-20, 20, by=1), 999)
          diffAccumColorsVec <- cm.colors(length(diffAccumBreaksVec)-1)
          agwPlotNiftyMatrix(recovery$accum
                             , mainLabel=agwGlue("Differential recovery of no. targeting drugs in top N (of ", length(targetingTestedDrugsVec), " possible drugs), vs. expected rate of ", format(chanceOfGettingPositiveRandomly, digits=2))
                             , imageColorsVec=diffAccumColorsVec
                             , imageBreaksVec=diffAccumBreaksVec)
          
          ## =======================================================================================
          screen(7)

          accumFracInTopMat <- numInTopMat
          for (topIdx in 1:ncol(numInTopMat)) {
               for (gene in nrow(numInTopMat):1) {
                    totalHittingDrugsInTop         <- sum(numInTopMat[nrow(numInTopMat):gene,topIdx])      #sum(numInTopMat[r, ncol(numInTopMat):c])
                    maxPossibleAcrossGenes <- ((nrow(numInTopMat)-gene)+1)*(min(numDrugsTargetingThisModule, topKVec[topIdx])) # how many top drugs COULD we have gotten?o
                    #maxPossibleAcrossGenes ## the total number of tested drugs that could even BE in the topk
                    accumFracInTopMat[gene, topIdx] <- totalHittingDrugsInTop/maxPossibleAcrossGenes #log(totalHittingDrugsInTop+0.1)
                    #accumFracInTopMat[r,c] <- sum(numInTopMat[nrow(numInTopMat):r,c])
               }
          }
          
          accumFracDivisions <- 10
          accumFracBreaksVec <- c(-999, 0, seq(1.0/accumFracDivisions,1.0,by=1.0/accumFracDivisions)) # -999 is for the "less than or equal to zero" range
          accumFracColorsVec <- agwGetHeatColors((length(accumFracBreaksVec)-1))
          agwPlotNiftyMatrix(accumFracInTopMat
                             , mainLabel=agwGlue("Fractional recovery of no. targeting drugs in top N (of ", length(targetingTestedDrugsVec), " possible drugs)")
                             , imageColorsVec=accumFracColorsVec
                             , imageBreaksVec=accumFracBreaksVec)
          
          close.screen(all=TRUE)
          agwFinishPlot()
          
          whichTopIdx <- 3
          ## plot(recovery$expected[,whichTopIdx], recovery$actual[,whichTopIdx], pch=' '
          ##      , asp=1.0
          ##      , main=agwGlue(length(recovery$expected[,whichTopIdx]), " : ", length(recovery$actual[,whichTopIdx]))
          ##      )
          ## abline(a=0, b=1, col="purple", lwd=2  )
          ## rect(xleft=-100, xright=0, ytop=999, ybottom=-999, col="#665566", border=NA)
          ## rect(xleft=-999, xright=999, ytop=0, ybottom=-999, col="#665566", border=NA)
          ## points(recovery$expected[,whichTopIdx], recovery$actual[,whichTopIdx], pch=PCH.SMALL.CIRCLE)
          ## lines(recovery$expected[,whichTopIdx], recovery$max[,whichTopIdx], pch=PCH.SMALL.CIRCLE, col="red", type='h')
          
          
          agwPreparePlot(directory=agwGlobalOutPath(agwGlue(k_FIGURE_3_FOLDER_NAME))
                         , file=agwGlue(agwAbridgeModuleNames(moduleName), "_recovery")
                         , res=100
                         , width=1024
                         , height=1024*2)
          close.screen(all=TRUE)
          split.screen(c(4,1))
          ## ==========================================

          miniGridColor      = "#00000033"
          bigGridColor       = "#66000066"
          underwaterColor    = "#00009966" # color for the "worse than random" region
          expectedPointPch   = '.'
          expectedPointColor = "gray"
          #actualColor = "orange"
          screen(1)
          
          indicatorGenesXlab = "No. indicator genes that were checked, ranked by PGI score"


          plot(c(recovery$expected[,whichTopIdx], max(recovery$actual[,whichTopIdx]), min(recovery$actual[,whichTopIdx])), pch=' ',
               , main=agwGlue(moduleName, ": Examining the top ", topKVec[whichTopIdx], " most-sensitizing drugs for each indicator gene.")
               , xlab="No. indicator genes used, ranked by PGI.\nLight green = max possible recovered drugs. Blue region = expected by chance."
               , ylab=agwGlue("No. *hitting* drugs recovered out of top ", topKVec[whichTopIdx], " most-sensitizing drugs."))
          
          polygon(x=c(0, length(recovery$expected[,whichTopIdx]), length(recovery$expected[,whichTopIdx])), y=c(0,max(recovery$expected[,whichTopIdx]),0), col=underwaterColor, border=NA) # "worse than random"
          points(rev(recovery$expected[,whichTopIdx]), pch=expectedPointPch, col=expectedPointColor)
          points(rev(recovery$actual[,whichTopIdx]), pch=PCH.DIAMOND,      col="black", bg=rev(sgiColorVec))
          points(rev(recovery$max[,whichTopIdx]), pch=PCH.X,            col="green") # max
          abline(h=seq(from=0, by=20, length.out=20), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=10, length.out=50), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=50, length.out=50), lwd=2, col=bigGridColor)
          ## ==========================================
          
          screen(2)
          plot(rev(sgiVec)
               , ylim=c(0, max(max(sgiVec)+0.1, 5))
               , xlab=indicatorGenesXlab
               , ylab="PGI Connectivity Score (-log10(p))"
               , pch=' ')
          agwFillRect(6000, "#FFEE9944")
          points(rev(sgiVec), pch=PCH.DIAMOND, col=rev(sgiColorVec), type='h')
          points(rev(sgiVec), pch=PCH.DIAMOND, col="black", bg=rev(sgiColorVec))
          abline(h=kSGI_CUTOFF_FIG3a, lwd=2, col=madeSGICutoffColor, lty="dashed")
          abline(h=seq(from=0, by=5, length.out=20), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=10, length.out=50), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=50, length.out=50), lwd=2, col=bigGridColor)
          ## ==========================================
          
          screen(3)
          plot2Thing <- rev(recovery$actual[,whichTopIdx]-recovery$expected[,whichTopIdx])
          plot(plot2Thing, pch=expectedPointPch, col=expectedPointColor
               , xlab=indicatorGenesXlab
               , ylab=agwGlue("No. *hitting* drugs - no. *expected*, out of top ", topKVec[whichTopIdx], " most-sensitizing drugs."))
          rect(xleft=-999, xright=999, ytop=0, ybottom=-999, col=underwaterColor, border=NA) # "worse than random" region
          points(rev(recovery$max[,whichTopIdx]-recovery$expected[,whichTopIdx]), pch=PCH.X, col="green") # max
          points(plot2Thing, pch=PCH.DIAMOND, col="black", bg=rev(sgiColorVec))
          abline(h=seq(from=0, by=10, length.out=20), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=10, length.out=50), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=50, length.out=50), lwd=2, col=bigGridColor)
          
          ## ==========================================
          screen(4)
          toPlot    <- rev(recovery$actual[,whichTopIdx]/recovery$max[,whichTopIdx])
          toPlotExp <- rev(recovery$expected[,whichTopIdx]/recovery$max[,whichTopIdx])

          toPlotModule    <- rev(recoveryInModule$actual[,whichTopIdx]/recoveryInModule$max[,whichTopIdx])
          toPlotExpModule <- rev(recoveryInModule$expected[,whichTopIdx]/recoveryInModule$max[,whichTopIdx])
          
          plot(c(toPlot, max(toPlotExp), min(toPlotExp)), pch=' '
               , ylim=c(0.0,1.0)
               , xlab=indicatorGenesXlab
               , ylab=agwGlue("Fraction of *hitting* drugs recovered, out of top ", topKVec[whichTopIdx], " most-sensitizing drugs."))
          #points(rev(recovery$max[,whichTopIdx]), pch=PCH.X, col="red") # max
          rect(xleft=-999, xright=999, ytop=max(toPlotExp), ybottom=-999, col=underwaterColor, border=NA) # "worse than random" region

          abline(h=mean(toPlotModule, na.rm=TRUE), lwd=4, col="red")
          #rect(xleft=-999, xright=999, ytop=max(toPlotExpModule), ybottom=-999, col=underwaterColor, border=NA) # "worse than random" region

          
          abline(h=seq(from=0, by=0.10, length.out=20), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=10, length.out=50), lwd=2, col=miniGridColor)
          abline(v=seq(from=0, by=50, length.out=50), lwd=2, col=bigGridColor)
          #rect(xleft=-100, xright=0, ytop=999, ybottom=-999, col="#665566", border=NA) # gray rect left of the Y-axis
          #rect(xleft=-999, xright=999, ytop=0, ybottom=-999, col="#665566", border=NA) # gray rect under the X-axis
          points(toPlotExp, pch=expectedPointPch, col=expectedPointColor)
          points(toPlot, pch=PCH.DIAMOND, col="black", bg=rev(sgiColorVec))

          #points(toPlotExpModule, pch=PCH.BOX, col="orange", cex=1)
          points(toPlotModule, pch=PCH.BOX, col="black", bg=rev(sgiColorVec), cex=1)
          
          ## ==========================================
          close.screen(all=TRUE)
          agwFinishPlot()
          #browser()
     }
}

agwFigure3a_main <- function() {
     print("=========================================")
     print("About to plot Figure 3A for all drugs...")
     print(agwGlue("Going to save in the Figure 3A directory named <", k_FIGURE_3_FOLDER_NAME, ">..."))
     agwFigure3aPlot()
}
