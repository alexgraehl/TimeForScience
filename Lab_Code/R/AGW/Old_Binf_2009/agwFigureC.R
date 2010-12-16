
##       source("/Users/alexgw/R/src/agwFigureC.R");        ## <-- To reload this file

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

OUTPUT_FIGURE_C_DIRECTORY = "~/R/FigureC_Drug_Clusters"

## Draws three figures in the OUTPUT_FIGURE_C_DIRECTORY directory.
## Clusters drugs that match certain names OR target names.
## IN OTHER WORDS, this does NOT show every single drug's information.
## It only shows drugs that either target a pathway in "targetPattern"
## or drugs whose names appear in "namePattern". Those are both search expressions,
## by the way.

agwFigure_C_ClusterDrugs <- function(drugOrfMatrix, targetPattern, namePattern, corType="spearman", outputLocation="~/R", outputPrefix="Plots") {
     agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
     agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")

     ## targetPattern: include drugs with an annotated TARGET that matches this perl regexp pattern. Example:   targetPattern="GoDDCR|DNA|RNA|mitochon"
     ## namePattern: include drugs with a NAME that matches this regexp.
     ## Figure out which drug names we are actually interested in showing on the heatmap,
     ## based on their annotated targets. We are selecting out ONLY drugs that target things that match
     ## lookForThisTargetPattern

     ## Add drugs that TARGET a specific regexp pattern (could be "DNA" or something like that)
     interestingDrugIndicesVec <- grep(targetPattern, GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)

     ## Add drugs whose NAMES match a specific regexp pattern
     additionalDrugIndicesVec <- grep(namePattern, colnames(drugOrfMatrix), ignore.case=TRUE, perl=TRUE)

     interestingDrugNamesVec <- union(rownames( GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ interestingDrugIndicesVec, ])
                                      , colnames(drugOrfMatrix[, additionalDrugIndicesVec ]))

     ## Below, commented out for now:
     ## This is a TEST, for debugging. You can see if this thing is working at all, by using this instead of the default gm calculations below
     ## If you switch the definition of "gm" to this stuff below,
     ## then you'll see how the plots work on sample data.
     #gm <- matrix(rnorm(100), nrow=10, byrow=T)
     #rownames(gm) <- paste("ROW", 1:10)  ;  colnames(gm) <- paste("COL", 1:10)
     #gm[12:30] = NA   ; gm[1:10] = 10:1   ;  gm[81:90] = c(10,9:1)
     #gm1 <- gm
     
     interestingDrugsNotInOurDataset <- setdiff(interestingDrugNamesVec, colnames(drugOrfMatrix))
     
     gm <- as.matrix(drugOrfMatrix[,
                                   intersect(colnames(drugOrfMatrix),  interestingDrugNamesVec) ] )
     
     agwPrint("Printing only the drugs which were annotated to target the pathways matching <", targetPattern, ">")
     agwPrint(ncol(gm), " drugs were found that met this requirement.")
     cat(length(interestingDrugsNotInOurDataset), "drugs were in the annotated-drugs set, but were not in our dataset. Here they are:\n")
     print(interestingDrugsNotInOurDataset)
     cat("\n")
     
     requireThisManyDataPointsPerDrug <- 3
     nColBefore <- ncol(gm)
     gm <- agwMatrixLinesWithEnoughData(gm, minN=requireThisManyDataPointsPerDrug) ## 3 or more data points required per row
     nColAfter <- ncol(gm)
     
     agwPrint("After requiring at least ", requireThisManyDataPointsPerDrug, " experimental data points per drug, we have ", nColAfter, " drugs remaining in the set to analyze. ", (nColBefore - nColAfter), " columns were removed.")

     ##gm[is.na(gm)] <- 0 ## dubious...
     
     USE_POLICY = "pairwise.complete.obs"
     drugCor <- cor(gm, method=corType, use=USE_POLICY)
#     orfCor  <- cor(t(gm), method=corType, use=USE_POLICY)

     agwPrint("Printing the distance matrix for the distance between the drugs...")
     drugDist <- as.dist(1-drugCor)
     agwPrint("Distance matrix has ", length(drugDist), " elements.")
#     orfDist  <- as.dist(1-orfCor)

     meanDrugDist <- mean(drugDist, na.rm=TRUE)
#     meanOrfDist  <- mean(orfDist , na.rm=TRUE)

     drugDist[is.na(drugDist)] <- meanDrugDist ## Hmmm... dubious!
#     orfDist[is.na(orfDist)]   <- meanOrfDist  ## Remove NAs for clustering

     agwPrint("Clustering the columns... (", length(drugDist), " data elements)")
     hColClust <- hclust( drugDist, method="complete", members=NULL)
     agwPrint("Done.")

#     agwPrint("Clustering the rows... (", length(orfDist), " data elements)")
     #hRowClust <- hclust( orfDist , method="complete", members=NULL)
     agwPrint("Done.")
     #par(mfrow = c(2, 2)); plot(hr, hang = 0.1); plot(hr, hang = -1)

     #browser()

#     m.top    = 2   ; m.left   = 2
#     m.bottom = 25  ; m.right  = 2
#     par(mar=c(m.bottom, m.left, m.top, m.right))

     z1 <- rep("", times=length(hColClust$labels))

     a <- grep("DNA Synth", GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)
     z1[ hColClust$labels %in% rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ a, ]) ] <- "### DNA-SYNTH-RELATED DRUG ### "

     a <- grep("Actin", GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)
     z1[ hColClust$labels %in% rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ a, ]) ] <- "===== Actin ===== "

     a <- grep("TOR Signaling", GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)
     z1[ hColClust$labels %in% rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ a, ]) ] <- "::::: TOR ::::: "

     a <- grep("Microtubule", GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)
     z1[ hColClust$labels %in% rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ a, ]) ] <- "<<<<< Microtubule >>>>> "
     
     a <- grep("Ergosterol", GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,3], ignore.case=TRUE, perl=TRUE)
     z1[ hColClust$labels %in% rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[ a, ]) ] <- "~~~~~ Ergosterol or Cell Wall ~~~~~ "

     z2 <- rep("", times=length(hColClust$labels))
     a <- grep("arabinof|hydroxy|fluoro|aminop|methotrex", hColClust$labels, ignore.case=TRUE, perl=TRUE)
     z2[ a ] <- ' ##### ANTIMETABOLITE ##'

     a <- grep("bleomycin|MMS|zocin|phleomycin|echino|cisplat|ellipti|actinomy|camptothe|neomycin|hygromyc", hColClust$labels, ignore.case=TRUE, perl=TRUE)
     z2[ a ] <- ' @@@@@ DIRECT DNA DAMAGE @@'
     
     hColClust$labels <- agwGlue(z1, hColClust$labels, z2)

     # =============================================
     agwPreparePlot(outputLocation, agwGlue(outputPrefix, "Dendrogram"), res=100, height=1500, width=500+(40*length(hColClust$labels)))
     plot(hColClust, lwd=6, col="black", xlab="", sub=agwGlue("Target pattern: <", targetPattern, ">. Name pattern: ", namePattern), main=agwGlue("Drugs clustered by degree of effect on growth of ORF deletion mutants. Correlation method: ", corType)
          )
     numClustersToShow  <- floor(length(hColClust$labels)/5)+1
     rect.hclust(hColClust, k=numClustersToShow, border="red") #rainbow(numClustersToShow))
     agwFinishPlot()
     # =============================================
     # =============================================
     agwPreparePlot(outputLocation, agwGlue(outputPrefix, "ColorDendrogram"), res=120, height=1500, width=500+(20 * length(hColClust$labels)))
     require(fpc)  ## Required for A2Rplot
     require(A2R)  ## Also required for A2Rplot <-- if you don't have this, you can actually just comment out the lines between here and "agwFinishPlot()", and you will get the other plots, just not the nice color dendrogram.
     m.top    = 2   ; m.left   = 2
     m.bottom = 40  ; m.right  = 2
     par(mar=c(m.bottom, m.left, m.top, m.right))
     A2Rplot(hColClust, k=numClustersToShow, boxes=TRUE
             , col.up="gray", col.down=rep(c("red","blue","orange","green"), numClustersToShow) #rainbow(numClustersToShow, v=0.8)
             , lwd.up=3, lwd.down=4
             , lty.up="solid", lty.down="solid"
             , type="rectangle", only.tree=TRUE, main="TITLE")
     mtext(hColClust$labels[hColClust$order], at=labels(hColClust$labels), side=1, outer=FALSE, las=2, line=2)
     agwFinishPlot()
     # =============================================
     # =============================================
     agwPreparePlot(outputLocation, agwGlue(outputPrefix, "Horizontal"), res=120, height=500+(20 * length(hColClust$labels)), width=2000)
     m.top    = 2  ; m.left   = 2
     m.bottom = 2  ; m.right  = 45
     par(mar=c(m.bottom, m.left, m.top, m.right))
     plot(as.dendrogram(hColClust), edgePar=list(col="blue", lwd=3), horiz=TRUE, type="rectangle")
     abline(v=c(0.5,1.0,1.5), col="orange", lty="dashed", lwd=1)
     agwFinishPlot()
     # =============================================
     
     agwPreparePlot(outputLocation, agwGlue(outputPrefix, "Heatmap"), res=120, height=4000, width=4000)

     agwPrint("Making a heatmap...")
     heatmap(t(gm[,]), Colv=NA  #Colv=as.dendrogram(hRowClust),
             , Rowv=as.dendrogram(hColClust), col=rainbow(20), scale="none"
             , margins=c(20,20))
     agwFinishPlot()
     agwPrint("Done! Made a heatmap...")

}


lookForTargets = "DNA Synth|HMG|Microtub|Ergosterol|Actin|Kinase Inhibitor|HDAC|Oxidation|TOR Signaling|Our_theoretical_tamox_target_check_for_cd"
lookForNames   = "benomyl|bleomycin|aminopter|methotrex|melphal"

## Use pearson correlation between drugs' vectors of growth inhibition in various ORF mutants
agwFigure_C_ClusterDrugs(drugOrfMatrix=GLOBAL.HOM.SENS.COMPENDIUM
                         , outputLocation=OUTPUT_FIGURE_C_DIRECTORY
                         , outputPrefix="Pearson_"
                         , targetPattern=lookForTargets
                         , namePattern=lookForNames
                         , corType="pearson")



## Use spearman correlation between drugs' vectors of growth inhibition in various ORF mutants
agwFigure_C_ClusterDrugs(drugOrfMatrix=GLOBAL.HOM.SENS.COMPENDIUM
                         , outputLocation=OUTPUT_FIGURE_C_DIRECTORY
                         , outputPrefix="Spearman_"
                         , targetPattern=lookForTargets
                         , namePattern=lookForNames
                         , corType="spearman")

