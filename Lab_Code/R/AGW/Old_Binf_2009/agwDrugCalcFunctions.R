## ================================================================

## source("/Users/alexgw/R/src/agwDrugCalcFunctions.R"); ## <-- To reload this file

source("/Users/alexgw/R/src/agwStandardFunctions.R");

DEFAULT.ROW.HEIGHT     = 14

HIST.DEFAULT.POINTSIZE = 18
hist.width             = 800
hist.default.height    = 800
hist.row.height        = DEFAULT.ROW.HEIGHT

boxplot.row.height = DEFAULT.ROW.HEIGHT   ## Height in pixels for each row
boxplot.width      = 1200  ## Width of the entire plot region
boxplot.pointsize  = 12   ## Font size
boxplot.m.top = 3; boxplot.m.left = 40; boxplot.m.bottom = 8; boxplot.m.right = 3;

dotchart.default.rowheight = DEFAULT.ROW.HEIGHT
dotchart.default.width     = 800
dotchart.default.pointsize = 12

DATA_ROOT_DIR               <- "~/R/Data"
PREDICTIONS_DIRNAME_IN_ROOT <- "Predictions"

SOURCE.HOM.COMPENDIUM.SENS.FILE <- agwGlue(DATA_ROOT_DIR, "/zz_hom_full.tab")
SOURCE.HET.COMPENDIUM.SENS.FILE <- agwGlue(DATA_ROOT_DIR, "/zz_het_full.tab")
KNOWN.TARGET.FILE               <- agwGlue(DATA_ROOT_DIR, "/zz_positives_no_orfs.tab")
KNOWN.TARGET.FILE.NO.DUPS       <- agwGlue(DATA_ROOT_DIR, "/zz_positives_no_orfs_or_redundancy.tab")
KNOWN.TARGET.FILE.ORFS          <- agwGlue(DATA_ROOT_DIR, "/zz_positives_ORF_set.tab")

gvGeneSLOverlap.FILENAME <- agwGlue(DATA_ROOT_DIR, "/gene_nonzero_sl_overlap.tab")
gvGeneSLPathway.FILENAME <- agwGlue(DATA_ROOT_DIR, "/module_nonzero_gene_sl_overlap.tab")

gvYeastGeneDescriptions.FILENAME <- agwGlue(DATA_ROOT_DIR, "/orf_descriptions.tab") ## Three columns: first is the ORF (sometimes missing), second is the gene name (always present), third is the description (a long text field)

GLOBAL.OUTPUT.PATH         = "~/R" ## <-- should *not* end in a '/'
ENRICH.LIMIT.SMALL.P.VALUE = 1e-17 ## 1e-17 seems to be the limit before enrich.R starts returning 0 as a P-value. (Obviously the p-value cannot actually be zero, so we need to clamp things that are less than this value to be this value.)




ADD.DEBUG.FAKE.DRUGS <- FALSE ## Set this to "True" to add a ridiculously good SGI-pathway-correlated set of data for testing purposes

foldDifferenceForSensitivityToCount <- 1 # <-- 1.0 means that one thing just needs to strictly be *more* sensitive in order to count
sensThresh <- log10(foldDifferenceForSensitivityToCount)
## Threshold for "ORF was differentially sensitive to drugs that targeted the pathway"
## In other words: a deletion orf counts as "differentially sensitive" for
## a specific pathway if the average growth for that ORF deletion when
## exposed to a drug is less than (foldDifferenceForSensitivityToCount) times less
## than the growth of the ORF compared to the average of all drugs.




## =================================================================

### Takes two vectors of (possibly) different length, and sticks them together (as columns)
### It's like cbind
##agwColBindNoRep <- function(c1, c2) {
##	
##	maxLen = max(length(c1), length(c2))
##	
##	newC1 = c(c1, rep(NA,maxLen-length(c1)))
##	newC2 = c(c2, rep(NA,maxLen-length(c2)))
##	
##	m <- matrix(data=c(newC1, newC2), byrow=FALSE, ncol=2)
##	return(m)
##}


## =================================================================

agwGetAnnotationFromOrf <- function(ORF) {
     if (!agwHasContent("gvYeastGeneDescriptions")) { agwGlobalLoad("gvYeastGeneDescriptions") }
     retVal <- (agwHashGet(hash=gvYeastGeneNameHash, key=ORF))
     if (is.null(retVal)) {
          return(list(gene=agwGlue(ORF," (No Name)")
                      , desc="(No description)"))
     } else {
          return(retVal)
     }
     # Note: returns a LIST, with the two sub-elements
     # "gene" (the gene name)
     # and "desc" (the long description)
}

## =================================================================

agwFindRanksForORFs <- function(compendium, byWhich) {
     ## Usage example: agwFindRanksForORFs("YPL233W", GLOBAL.HOM.SENS.COMPENDIUM)

     ## Note! This is a very slow function, because it has to sort a lot of data!
     ## Orf: a textual ORF name like "YPL233W"
     ## Compendium: a 2-D data frame in the style of GLOBAL.HOM.SENS.COMPENDIUM
     ##  (it is drugs in rows vs. orfs in columns,
     ## with each cell being a sensitivity score number)

     ## Note that this is the ORFs with their sensitivies replaced by rank.
     ## If you want to list the DRUGS in order of rank, then this is NOT the
     ## function you want!

     ## byWhich is "rank by row, or by column?" Column = by DRUGS, Row = by ORFs
     print("Running agwFindRanksForORFs... this will take a minute!")

     ## Sort from high to low:
     sortFromHighToLow <- TRUE
     rankedMatrix <- apply(compendium, byWhich, function(oneLineVec) {
          cat("~") ;
          z <- rank((-1 * oneLineVec), na.last=TRUE, ties.method="average"); ## <-- note, we NEGATE a here, because we want to rank from high values (good ranks: 1,2,3,...) to low/negative values (bad ranks: 99,100,101,...). Also: NAs are NOT counted as equal for ranks!!!
          numNonNA <- sum(!is.na(oneLineVec))
          z[z > numNonNA] <- NA ## NA input data gets NA ranks
          ##print(numNonNA);
          return(z)
     } )
     
     print("\nDone running agwFindRanksForORFs...")

     return(rankedMatrix)
     ## rankedMatrix is like the input, except the values have been all
     ## replaced with their ranks instead!
     ## An input matrix with values  888, 11, 0, NA, -999, 199999999
     ## would be converted to          2,  3, 4,  5,    6, 1
     ## Ties are given the AVERAGE of the tying elements.
     
     ## This function seems to actually work. To verify, check out:
     ## drug = 4 ; orf = 2:10; GLOBAL.HOM.SENS.COMPENDIUM[orf,drug] ; gTest[orf,drug] ; colnames(GLOBAL.HOM.SENS.COMPENDIUM)[drug] ; colnames(gTest)[drug] ; rownames(GLOBAL.HOM.SENS.COMPENDIUM)[orf] ; rownames(gTest)[orf]
     
}


## =================================================================
agwGlobalFileLoad <- function(varname) {
     return(agwGlobalLoad(varname))
}

agwGlobalLoad <- function(varname) {
     ## You pass in a STRING, indicating the variable that you want
     ## to have loaded / populated. If this global variable does
     ## not already exist, it will be added to the global namespace
     ## (using the operator "<<-" )
     if (agwHasContent(varname)) {
          ## Variable is *already* defined, AND has some content (i.e. is not an empty list)
          cat(agwGlue(">>>> agwGlobalLoad: It isn\'t necessary to reload <", varname, ">, which was loaded earlier.\n"))
          return(FALSE) ## did not load new data...
     } else {
          cat(agwGlue(">>>> agwGlobalLoad: About to load data for the variable <",varname,"> into the global namespace...\n"))
     }
     
     if (varname == "GLOBAL.HOM.SENS.COMPENDIUM") {
          GLOBAL.HOM.SENS.COMPENDIUM <<- agwReadFileIntoDataFrame(SOURCE.HOM.COMPENDIUM.SENS.FILE)

     } else if (varname == "GLOBAL.HET.SENS.COMPENDIUM") {
          GLOBAL.HET.SENS.COMPENDIUM <<- agwReadFileIntoDataFrame(SOURCE.HET.COMPENDIUM.SENS.FILE)
          
     } else if (varname == "GLOBAL.POSITIVES.COMPLETE.DATAFRAME") {
          GLOBAL.POSITIVES.COMPLETE.DATAFRAME <<- agwReadFileIntoDataFrame(KNOWN.TARGET.FILE)
          
     } else if (varname == "KNOWN.TARGETS.WITH.ORFS") {
          KNOWN.TARGETS.WITH.ORFS <<- agwReadFileIntoListOfLists(KNOWN.TARGET.FILE.ORFS, row.name.column=1)

     } else if (varname == "gv.HOM.RANKS.BY.DRUG") {
          ## NOTE! Not the same as the "by orf" variant
          agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")          
          gv.HOM.RANKS.BY.DRUG <<- agwFindRanksForORFs(GLOBAL.HOM.SENS.COMPENDIUM, BY.COL) ## For each drug, the ORF sensitivities are replaced by the RANK (1 = most sensitive)
          
     } else if (varname == "gv.HOM.RANKS.BY.ORF") {
          ## Ranks of each DRUG for each ORF. Not the same as gv.HOM.RANKS
          agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
          gv.HOM.RANKS.BY.ORF <<- agwFindRanksForORFs(GLOBAL.HOM.SENS.COMPENDIUM, BY.ROW) ## For each drug, the ORF sensitivities are replaced by the RANK (1 = most sensitive)

     } else if (varname == "GLOBAL.POSITIVES.NO.DUPS.DATAFRAME") {
          GLOBAL.POSITIVES.NO.DUPS.DATAFRAME <<- agwReadFileIntoDataFrame(KNOWN.TARGET.FILE.NO.DUPS)

     } else if (varname == "gvGeneSLPathway") {
          gvGeneSLPathway <<- agwReadFileIntoDataFrame(gvGeneSLPathway.FILENAME, allowRagged=TRUE, row.names=NULL)
          
     } else if (varname == "gvGeneSLOverlap" || varname == "gvGeneSLOverlap.HASH") { ## Amount of overlap of a gene's SL neighbors with our modules
          gvGeneSLOverlap <<- agwReadFileIntoDataFrame(gvGeneSLOverlap.FILENAME, allowRagged=TRUE, row.names=NULL)
          ## Creates a new "environment" named "gvGeneSLOverlap.HASH" ("SL hash"). An environment is kind of like a hash,
          ## in that you can access things in it in constant time, and store keys/values.               
          ## The only problem is that it's fairly slow to populate. (~1 minute for 300,000 entries)
          ##slh <- list()
          gvGeneSLOverlap.HASH <<- agwNewHash()
          totalRows = nrow(gvGeneSLOverlap)
          for (i in 1:totalRows) {
               gene      = gvGeneSLOverlap[["Deletion_Gene"]][i]
               mod       = gvGeneSLOverlap[["Module"]][i]
               agwHashPut(hash=gvGeneSLOverlap.HASH, key=paste(gene,mod,sep="~"), value=i) ## <-- save the index of the value in the original array
               if (i %% 10000 == 0) {
                    print(paste("Populating the overlap score hash... done with ", i, "out of", totalRows))
               }
          }
          print(paste("Done reading",length(gvGeneSLOverlap.HASH),"elements into gvGeneSLOverlap.HASH."))
          
     } else if (varname == "gvYeastGeneDescriptions" || varname == "gvYeastGeneNameHash") {

          gvYeastGeneDescriptions <<- agwReadFileIntoDataFrame(gvYeastGeneDescriptions.FILENAME, allowRagged=TRUE, row.names=NULL)
          gvYeastGeneNameHash     <<- agwNewHash()
          for (i in 1:nrow(gvYeastGeneDescriptions)) {
               agwHashPut(hash=gvYeastGeneNameHash
                          , key=gvYeastGeneDescriptions[,1][i]
                          , value=list(gene=gvYeastGeneDescriptions[,2][i], desc=gvYeastGeneDescriptions[,3][i]) )
          }
          
     } else if (varname == "gvMemberStatsList" || varname == "gvNeighborStatsList" || varname == "gvNeighborStatsHash" || varname == "gvModVec") {
          agwSetUpDiffSensVsCentralityScoreGlobalVariables() ## sort of hack-ish... has to calculate a bunch of stuff for each drug. This is defined in agwDrugCalcFunctions.R
     } else if (varname == "gvCompleteCollection") {
          agwSrcAndRun("/Users/alexgw/R/src/agwFigureA_DetailPlots.R")
          # this should load things into gvCompleteCollection...
          
     } else if (varname == "gvDrugTargetMapping") {
          ## Creates a "tarlist" target list. This is just a list of lists telling us
          ## which drugs were listed as targeting which modules.
          ## The order is:  tarlist[[MODULE_NAME]][target_name]
          
          agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME") ## <-- required
          gvDrugTargetMapping <<- list() ## <-- global namespace
          print("gvDrugTargetMapping: Populating the drug target table...")
          drugsWithTarget  <- rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME)
          targetedModules  <- GLOBAL.POSITIVES.COMPLETE.DATAFRAME[["TARGET_MODULE"]]
          tarModLen <- length(targetedModules)
          for (i in 1:tarModLen) {
               tar <- drugsWithTarget[i]
               mo  <- targetedModules[i]
               ## For each module ("mo"), figure out which drugs are listed to target that pathway
               ## Sets a location in the list to "1" for every pair where the drug is listed to target the pathway.
               if (is.null(gvDrugTargetMapping[[mo]][tar]) || is.na(gvDrugTargetMapping[[mo]][tar])) {
                    gvDrugTargetMapping[[mo]][tar] <<- 1
               }
          }

     } else if (varname == "gvPathwayMembershipHash") {
          ## Note that pathway membership is shown ONLY for drugs with known targets.
          ## This is not a list of *all* known targets.
          gvPathwayMembershipHash <<- agwNewHash()
          pathwayIdx  <- 4   ;    orfStartIdx <- 7 ## <-- column indices in the data structure
          agwGlobalLoad("KNOWN.TARGETS.WITH.ORFS")
          for (i in 1:length(KNOWN.TARGETS.WITH.ORFS)) {
               agwHashPut(key=KNOWN.TARGETS.WITH.ORFS[[i]][pathwayIdx]
                          , value=KNOWN.TARGETS.WITH.ORFS[[i]][orfStartIdx:length(KNOWN.TARGETS.WITH.ORFS[[i]])]
                          , hash=gvPathwayMembershipHash)
          }
          
     } else {
          stop(paste("Sorry, the agwGlobalLoad function (in agwDrugCalcFunctions.R) doesn't know about this variable name:", varname))
     }

     cat(agwGlue(">>>> Finished loading <",varname,"> into the global namespace.\n"))
     return(TRUE) ## loaded new data...
}
## =================================================================


## =================================================================
agwHandleDifferentialSensitivity <- function() {

     ## Calculates the "differential sensitivity" of each deletion-ORF mutant.
     ## Differential sensitivity is basically "how much more sensitive was this ORF mutant
     ## to drugs of a certain drug group / drugs with certain annotated targets, as compared
     ## to its sensitivity to other drugs"
     ## The actual mathematical definition is somewhat more complicated.

     agwGlobalLoad("gvGeneSLPathway") ## <-- this file is already in order from best-to-worst scores for each gene
     ##hist(gvGeneSLPathway[["SL_Overlap"]])	
     
     maxIndicatorGenes = 5 ## <-- How many top "indicator genes" do we want to look at for each pathway?
     
     ## Note: this is fairly slow (takes about 30 seconds on the Go set I'm using
     delVec <- NULL ## delVec is a vector with (vectors of gene modules) as elements. 
     ## In other words:  delVec[["some module"]] will return c(42349,2421,2425), where the numbers in that vector
     ## are the indices of the ORFs in the huge "gvGeneSLPathway" list that are potential indicator genes (i.e., the ORFs
     ## that have high synthetic lethal overlap with the module in question)
     
     if (!agwHasContent("delVec")) {
          delVec <- vector("list") ## <-- assign to the global-scope delVec that gets passed in
          numRows <- nrow(gvGeneSLPathway)
          DELETION_GENE_INDEX = 1
          MODULE_NAME_INDEX   = 2
          for (i in 1:numRows) {
               thisGene = gvGeneSLPathway[[DELETION_GENE_INDEX]]
               thisMod = gvGeneSLPathway[[MODULE_NAME_INDEX]][i]
               if (is.null(thisMod) || (nchar(thisMod) == 0)) { print("Skipping empty module") ; next ; } ## skip this iteration of the "for" loop
               if (is.null(delVec[[thisMod]])) { delVec[[thisMod]] = vector() }
               
               vectorSize <- length(delVec[[thisMod]])
               if (vectorSize < maxIndicatorGenes) {
                    delVec[[thisMod]] <- append(delVec[[thisMod]], i) ## store the current index
               }
               if (i %% 25000 == 0) { print(paste("Populating delVec: Done with",i,"of", numRows,"runs.")) }
          }	
          print(paste("Status report: read in a total of this many pathways:",length(delVec)))
     }
     
     agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")
     agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
     agwGlobalLoad("GLOBAL.POSITIVES.NO.DUPS.DATAFRAME")
     agwGlobalLoad("gvDrugTargetMapping")
     
     kORF.SEP.CHARS = " :: "
     
     agwGlobalLoad("gvYeastGeneDescriptions")
     kDESC.GENE.COL = 2 ## column 2 in the gvYeastGeneDescriptions file has gene names
     kDESC.INFO.COL = 3 ## column 3 in the gvYeastGeneDescriptions file has the long descriptions
     
     drugsThatWereTested = colnames(GLOBAL.HOM.SENS.COMPENDIUM)
     tTestInfoList = list() ## will save all the t-test info
     for (i in 1:length(gvDrugTargetMapping)) {
          mo <- names(gvDrugTargetMapping[i]) ## <-- a single module name

          stopifnot(length(mo) == 1) ## <-- sanity check to make sure we just got a single name back
          drugsThatTargetMo            <- names(gvDrugTargetMapping[[i]]) ## <-- a vector of *all* the drugs that are listed to target "mo"
          drugsThatDoNotTargetMo       <- setdiff(rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME), drugsThatTargetMo)
          testedDrugsThatTargetMo      <- intersect(drugsThatTargetMo, drugsThatWereTested)
          testedDrugsThatDoNotTargetMo <- intersect(drugsThatDoNotTargetMo, drugsThatWereTested)
          
          ##print("Tested drugs that target the module:") ; print(testedDrugsThatTargetMo)
          ##print(agwGlue("The known target of the current drug is <", mo,">."))
          ##print("This pathway is covered the most by the SL neighborhood of these ORFs:")
          
          for (idx in delVec[[mo]]) {
               ## The IDX runs through the various indices for ORFs that had significant overlap with the pathway in question
               theORF <- gvGeneSLPathway[["Deletion_Gene"]][idx]
               
               sensTargeted    <- GLOBAL.HOM.SENS.COMPENDIUM[theORF, testedDrugsThatTargetMo]
               sensNotTargeted <- GLOBAL.HOM.SENS.COMPENDIUM[theORF, testedDrugsThatDoNotTargetMo]
               
               theList <- list(as.vector(sensTargeted,mode="double"), as.vector(sensNotTargeted,mode="double"))
               numHitsNonNA  = length(na.omit(theList[[1]]))
               numOtherNonNA = length(na.omit(theList[[2]]))
               
               if (numHitsNonNA >= 2 && numOtherNonNA >= 2) { ## T test requires 2 obs. for both sets
                    ttt <- t.test(sensTargeted, sensNotTargeted,alternative="greater")
                    theRealGeneName <- gvYeastGeneDescriptions[theORF, kDESC.GENE.COL]
                    theDescription  <- gvYeastGeneDescriptions[theORF, kDESC.INFO.COL]
                    theIndex = length(tTestInfoList)+1
                    tTestInfoList[[theIndex]] <- ttt
                    tTestInfoList[[theIndex]]$`fullname` <- agwGlue( agwPerlGSubi("OUR_KNOWN_TARGET:","",mo), kORF.SEP.CHARS, theORF, kORF.SEP.CHARS, theRealGeneName," (p=",round(ttt$p.value,digits=3),", t=",round(ttt$statistic,digits=2),")")
                    ##print(agwGlue("Module ##", i, " <", mo, "> and ORF <", theORF, ">, ", ttt$method, "(alt: ", ttt$alternative, ")", " has t-stat ", round(ttt$statistic,digits=2), " and p-value ", round(ttt$p.value,digits=2)))
                    
                    ##print("SensTar") ; print(sensTargeted)
                    ##print("SensNot") ; print(sensNotTargeted)
                    
                    additionalTitle = agwGlue(theRealGeneName," (", theORF,"), p=",round(ttt$p.value,digits=5),", t=",round(ttt$statistic,digits=2))
                    
                    FILE = agwGlobalOutPath(agwGlue("T_TestDoubleBar/DoubleBar_", agwAbridgeModuleNames(mo),"_",theORF ))
                    agwHistCanvas(FILE, agwDoubleBox, theORF, mo, additionalTitle, sensTargeted, sensNotTargeted, testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo)
               } ## <-- end of "if the t-test is do-able"
          } ## end of looping through all the modules with known target/indicators ORFs AND not too many targets	
     }
     
     pDisplayMinCutoff = 0.05
     pDisplayMaxCutoffInLogSpace = 10  ## don't bother distinguishing between 10e-10 and further
     
     vp <- vector()
     vname <- vector()
     for (i in 1:length(tTestInfoList)) {
          ##vp <- append(vp, tTestInfoList[[i]]$p.value)
          ##vp <- append(vp, tTestInfoList[[i]]$statistic)
          loggedP <- -(log10(tTestInfoList[[i]]$p.value))
          clippedP <- min(loggedP, pDisplayMaxCutoffInLogSpace)
          vp <- append(vp, clippedP)
          vname <- append(vname, agwGlue(tTestInfoList[[i]]$`fullname`, ", -logP=", round(loggedP,digits=1)))
     }
     
     ##vname <- sort(vname) ## <-- this breaks things for some reason!
     ##vp <- vp[order(vname)]
     
     
     vnameModuleOnly <- agwPerlGSubi(agwGlue(kORF.SEP.CHARS,".*"),"", vname)  ## Used to compute the color "hash" for each module
     fullColorVec <- agwColorsFromStrings(vnameModuleOnly, saturation=1, value=0.75)  ## Color for points
     bgColorVec <- agwColorsFromStrings(vnameModuleOnly, saturation=0.35, value=1)    ## Color for background
     rm(vnameModuleOnly)
     

     plotMaxX = max(vp)+10 ## Max X value is the maximum p-value plus 10 (arbitrarily-chosen), just so that the labels fit on the plot

     FILE = agwGlobalOutPath("T_TestDoubleBar_Summary.png")
     if (!is.null(FILE)) { png( FILE, pointsize=HIST.DEFAULT.POINTSIZE, width=STANDARD.PLOT.WIDTH*3, height=STANDARD.PLOT.HEIGHT*2) }


     plot(x=vp, y=c(1:length(vp)), yaxt="n", xlim=c(0,plotMaxX), pch='.'
          ,xlab="-log(10) T-test p-value (0 = worst) for separation (> only) of the DeletionORF mutant's sensitivity to\nmodule-targeting drugs vs non-module-targeting drugs"
          ,ylab=agwGlue("Each indicator gene for a module (up to ", maxIndicatorGenes, " per module) (N = ",length(vp),")")
          ,main=agwGlue(
           "Differential growth of indicator DeletionORFs (chosen for pathway SGI coverage).\nQ: Are indicators more sensitive to drugs that target the pathway?")
          ,sub=agwGlue("Gray area is p > ", pDisplayMinCutoff," (uncorrected). X values capped at ", pDisplayMaxCutoffInLogSpace,"."))

     charExpansionSizer = 1.0
     rectXMin = 0  ;  rectXMax = pDisplayMaxCutoffInLogSpace ##plotMaxX
     rect(rectXMin, 0.5+seq(0,length(vp)-1), rectXMax, 0.5+seq(1,length(vp))+0.25, col=bgColorVec,lwd=0)
     rect(rectXMin, 0.5, -log10(pDisplayMinCutoff), (0.5+length(vp)+0.25), col="gray",lwd=0)
     text(x=vp, y=c(1:length(vp)), labels=vname, cex=charExpansionSizer, offset=charExpansionSizer, pos=POS.RIGHT)
     points(x=vp, y=c(1:length(vp)), col="black", bg=fullColorVec,pch=PCH.DIAMOND)

     if (!is.null(FILE)) { dev.off() }


}



## =================================================================

agwPredBoxPlot <- function(inMatrix, outFilename=NULL, numModulesTested=NULL) {

     pred.multiple.testing = ifelse(is.null(numModulesTested), nrow(inMatrix), numModulesTested) ## how many modules were tested as possible "hits" per drug?

     confidence.0.05  = (-log10(0.05 / pred.multiple.testing))
     confidence.0.01  = (-log10(0.01 / pred.multiple.testing))
     confidence.0.001 = (-log10(0.001 /pred.multiple.testing))

     mtext.line = 6 ## <-- for formatting the legend that goes by the x axis
     pred.x.label = agwGlue(
     "-log10(p-val) (Uncorrected for multiple testing).\n(0 = no prediction, right = more confident prediction)\n",
     "Bonferroni-corrected cutoffs (for ", pred.multiple.testing, " tests) shown at p=0.05, 0.01, and 0.001.\n",
     "Table dimensions: ",nrow(inMatrix)," rows by ",ncol(inMatrix), " columns.");
     pred.pch = PCH.DIAMOND
     pred.outlier.color = "red"

     if (!is.null(outFilename)) {
          png( outFilename,pointsize=boxplot.pointsize, width=boxplot.width, height=(100 + boxplot.row.height * ncol(inMatrix)) )
     }
     par(mar=c(boxplot.m.bottom, boxplot.m.left, boxplot.m.top, boxplot.m.right))
     boxplot(inMatrix, horizontal=TRUE, xlab='', main=paste("Predictions"), pch=pred.pch, bg=pred.outlier.color, col=agwColorsFromStrings(colnames(inMatrix)), axes=TRUE,las=LAS.HORIZONTAL.TEXT, notch=FALSE, border=c("black"))

     agwLine(confidence.0.05,'v') ; agwLine(confidence.0.01,'v') ; agwLine(confidence.0.001,'v')

     mtext(pred.x.label, side=POS.BOTTOM, outer=FALSE,line=mtext.line)
     if (!is.null(outFilename)) {
          dev.off()
     }

}
## =================================================================

agwGlobalOutPath <- function(...) {
     return (paste(GLOBAL.OUTPUT.PATH,"/", ..., sep=''));
}

## =================================================================

agwHistCanvas <- function(out, func, ...) {
     if (!is.null(out)) {
          png(filename=agwGlue(out,".png"), pointsize=HIST.DEFAULT.POINTSIZE, width=hist.width, height=hist.width)	}
     func(...) ;
     if (!is.null(out)) { dev.off(); }
}

## =================================================================

agwDoubleBox <- function(theORF, modname, additionalTitle, sensTargeted, sensNotTargeted, testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo) {
     ## This is a double box plot graph. One box goes on the left, and one goes on the right.
     ## The one on the left is for the drugs that claim to target a particular pathway, and the
     ## one on the right is for the drugs that claim to target some OTHER pathway.
     ## We are measuring the difference between them.

     ## must be called from WITHIN agwHistCanvas
     ## Example:
     ## FILE = agwGlobalOutPath(agwGlue("DoubleBar_", agwAbridgeModuleNames(mo)))
     ## agwHistCanvas(FILE, agwDoubleBox, mo, sensTargeted, sensNotTargeted, testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo)

     abridgedMo <- agwAbridgeModuleNames(modname) ## "abridged module name"

     colorSet = c("red","gray")

     ##print("Sens1") ; print(sensTargeted) ; 	print(sensNotTargeted) ;	print("Sens2")
     theList <- list(as.vector(sensTargeted,mode="double"), as.vector(sensNotTargeted,mode="double"))

     names(theList) <- c(paste(length(testedDrugsThatTargetMo), "drugs that target this module"), paste(length(testedDrugsThatDoNotTargetMo), "drugs with other targets"))
     ##print("Length of list:") ; print(length(theList))
     ##print(length(na.omit(theList))) ; print("The list:")
     ##print(theList) ; print("Done with list")

     numHitsNonNA  = length(na.omit(theList[[1]]))
     numOtherNonNA = length(na.omit(theList[[2]]))

     if ((numHitsNonNA + numOtherNonNA) > 0) {
          ymin = min(theList[[1]], theList[[2]],na.rm=TRUE)
          ymax = max(theList[[1]], theList[[2]],na.rm=TRUE)
          if (ymin > (-ymax)) { ymin = -ymax } ## "centralizes" the ylimits so that they are the same...
          if (ymax < (-ymin)) { ymax = -ymin } ## ...distance from zero

          boxplot(theList, col=colorSet, bg=colorSet, pch=PCH.DIAMOND, xlab="", ylab="log10(Growth_control/Growth_expy): higher = drug is more deadly", main=agwGlue(abridgedMo," :: ",theORF,"\n",additionalTitle), ylim=c(ymin,ymax))
          points(jitter(rep(1.5,length(theList[[1]])),0.25), theList[[1]], pch=4)
          points(jitter(rep(2.5,length(theList[[2]])),0.25), theList[[2]], pch=4)
          agwLine(0,'h')
     }
}


## =================================================================

agwHistogram <- function(inData, breaks="Sturges", out=NULL, xlab=NULL, title=NULL) {
     if (is.null(title) && !is.null(xlab)) {
          title <- paste("Distribution of", xlab);
     }

     if (!is.null(out)) {
          png( out, pointsize=HIST.DEFAULT.POINTSIZE, width=hist.width, height=(100 + hist.row.height * ncol(inData)) )
     }
     ##par(mar=c(boxplot.m.bottom, boxplot.m.left, boxplot.m.top, boxplot.m.right))
     hist(inData, breaks=breaks, col="lightblue", labels=TRUE, main=title, xlab=xlab, freq=TRUE)
     rug(inData, col="blue")

     ##agwLine(confidence.0.05,'v') ; agwLine(confidence.0.01,'v') ; agwLine(confidence.0.001,'v')
     ##mtext(pred.x.label, side=POS.BOTTOM, outer=FALSE,line=mtext.line)
     if (!is.null(out)) { dev.off() }
}
## =================================================================

## =================================================================

agwAbridgeModuleNames <- function(drugNames) {
     abridgedLabels = agwPerlGSubi("SL_Neighbors_of:","SL:", drugNames)
     abridgedLabels = agwPerlGSubi("OUR_KNOWN_TARGET:","KNOWN:", abridgedLabels)
     abridgedLabels = agwPerlGSubi("DRUGBANK_MAPPED_FROM_HUMAN:","DB_FROM_HUMAN:", abridgedLabels)
     abridgedLabels = agwPerlGSubi("AGW_Target_from_Imming","AGW_Imming", abridgedLabels)
     return(abridgedLabels)
}



agwSlidingWindowSmoother <- function(x,y, windowSize=1, numSteps=20) {
     ### Input: vectors x and y (typically these are plottable coordinates)
     ### windowSize--the size of the sliding window to average over (note that
     ### this is the TOTAL window size (the "diameter" of the window). It is not
     ### the "radius" of points that are included in a center point.)
     
     ### numSteps: number of steps from the beginning to the end

     ### Output: a list with:
     ### a new x coordinate (at the same location as the input x)
     ### a new y coordinate (the "smoothed" y value)
     ### an "n" (the count of items that went into calculating the y value)

     ### A manual version of "supsmu" or "loess.smooth"
     ### Note that it just does a sliding window, so sometimes the results
     ### are wonky. There is no "smoothing" function, it's just a "cliff"!

     ## Example usage:
     ### xpts <- c(1,2,4,5,6,4,5,9,2,4,10,8,12,4,3,5,4)
     ### ypts <- c(9,1,2,3,1,4,3,5,4,3,13,4,55,1,2,3,4)
     ### smoothedLineCoords <- agwSlidingWindowSmoother(xpts, ypts, windowSize=2, numSteps=20)
     ### plot(smoothedLineCoords$x, smoothedLineCoords$y, type='l')

     stopifnot(numSteps >= 2) ; stopifnot(windowSize >= 0) ; stopifnot(length(x) == length(y))
     windowStart  <- agwMin(x) ;  windowEnd <- agwMax(x)
     stepSize     <- (windowEnd - windowStart)/numSteps
     smoothedY    <- list()
     halfWindow   <- windowSize/2.0
     for (centerLoc in seq(from=windowStart, to=windowEnd, by=stepSize)) {
          inRangeVec <- (x >= (centerLoc-halfWindow) & x < (centerLoc+halfWindow))
          theCount <- sum(inRangeVec)

          if (theCount > 0) {
               theAvg <- mean(y[inRangeVec], na.rm=TRUE)
               ##agwPrint("Step ", round(stepSize,2), ", range ", round(xLeft,1), "---", round(xLeft+windowSize,1), ", adding the count ", thisCount, " (", thisSum,")--",round(thisSum/thisCount,2))
               idx = length(smoothedY) + 1
               smoothedY[[idx]] <- list(     avg = theAvg
                                        , center = centerLoc
                                        ,  count = theCount)
          }
     } ## <-- End of FOR loop
     ## Note: returns a list with an "x" and "y" component--easy for plotting with "lines"
     return(list(  x=sapply(smoothedY,"[[","center")
                 , y=sapply(smoothedY,"[[","avg")
                 , n=sapply(smoothedY,"[[","count")))
     ## This "rect" stuff below shows the sliding windows for the manually-generated smoothing function
     ##   rect(xleft=seq(from=windowStart-halfWindow, to=windowEnd-halfWindow, by=stepSize) ## Show the window size
     ##        , ybottom=seq(from=0,to=1,by=1/numSteps)
     ##        , xright=seq(from=windowStart+halfWindow, to=windowEnd+halfWindow, by=stepSize)
     ##        , ytop=seq(from=0,to=1,by=1/numSteps)+0.01
     ##        , col="yellow")
     ## abline(v=seq(from=windowStart, to=windowEnd, by=stepSize),col="green")
}


agwSlidingWindowSmootherCumulative <- function(xVec, yVec) {
     ### Input: vectors x and y (plot coordinates)
     ### numSteps: number of steps from the beginning to the end
     ### Note: this is a CUMULATIVE-style version of the agwSlidingWindowSmoother
     ### i.e., it starts from the highest, then moves its way down, averaging over all the
     ### previous data points. It has not been optimized for speed at all.

     ### Output: a list with:
     ### a new x coordinate (at the same location as the input x)
     ### a new y coordinate (the "smoothed" y value)
     ### an "n" (the count of items that went into calculating the y value)

     stopifnot(length(xVec) == length(yVec))
     sortedXVec    <- sort(xVec)
     ySortedByXVec <- yVec[ order(xVec) ]
     smoothedList  <- list()
     theSum   <- 0
     theCount <- 0
     for (i in length(sortedXVec):1) {
          theSum   <- (theSum + ySortedByXVec[i])
          theCount <- (theCount + 1)
          idx = length(smoothedList) + 1
          smoothedList[[idx]] <- list(  x = sortedXVec[i]
                                      , y = (theSum/theCount) ## <-- mean value
                                      , n = theCount)
     }
     return(list(  x=sapply(smoothedList,"[[","x")
                 , y=sapply(smoothedList,"[[","y")
                 , n=sapply(smoothedList,"[[","n")))
     ## Note: returns a list with an "x" and "y" component--easy for plotting with "lines"
}


agwFigureOutGeneInfo <- function(orf, indexInSLCollection, drugsThatTarget, drugsThatDoNot) {
     sensTargetedVec    <- as.vector(GLOBAL.HOM.SENS.COMPENDIUM[orf, drugsThatTarget], mode="double")
     sensNotTargetedVec <- as.vector(GLOBAL.HOM.SENS.COMPENDIUM[orf, drugsThatDoNot] , mode="double")
     rankTargetedVec    <- NA
     rankNotTargetedVec <- NA
     
     rowIndexForThisOrf <- match(orf, rownames(gv.HOM.RANKS.BY.DRUG)) # finds the proper row with the data for this particular ORF in it. I believe this is always a SINGLE value, but match can hypothetically return more than one.
     stopifnot(length(orf) == 1) # <-- if there is more than one row with the same ORF value in it, then things are completely messed up! This data is in a table, and there should only be one entry per orf!
     if (!is.na(rowIndexForThisOrf)) {
          rankTargetedVec    <- gv.HOM.RANKS.BY.DRUG[rowIndexForThisOrf, drugsThatTarget] ## <-- ranks! Not values!
          rankNotTargetedVec <- gv.HOM.RANKS.BY.DRUG[rowIndexForThisOrf,  drugsThatDoNot]  ## <-- ranks! Not values!
     }
     
     ttt <- NULL
     if (length(na.omit(sensTargetedVec)) > 0
         && length(na.omit(sensNotTargetedVec) > 0)) {
          theRealGeneName <- agwGetAnnotationFromOrf(orf)$gene
          
          overlapAmt <- NA
          if (!is.null(indexInSLCollection)) {
               ## Note: the line below works without [[1]]
               overlapAmt <- gvGeneSLPathway[indexInSLCollection,]["SL_Overlap"][[1]]
          } else {
               overlapAmt <- 0.0 ## No SGI prediction overlap found for this gene with the target module in question
          }

          theSensMedAbsDev <- mad( c(sensTargetedVec, sensNotTargetedVec)
                                  , constant=1, na.rm=TRUE) # mad is "median absolute deviation" (http://en.wikipedia.org/wiki/Median_absolute_deviation)

          theRankMedAbsDev <- mad( c(rankTargetedVec, rankNotTargetedVec)
                                  ,  constant=1, na.rm=TRUE) # mad is "median absolute deviation" (http://en.wikipedia.org/wiki/Median_absolute_deviation)

          #browser()
          ttt <- list(    rankMedIN = median(rankTargetedVec   , na.rm=TRUE) # <-- ranks by drug sensitivity (I think 1=most sens?)! not values
                      ,  rankMedOUT = median(rankNotTargetedVec, na.rm=TRUE) # <-- ranks, not sensitivities!
                      ,  rankMeanIN =   mean(rankTargetedVec   , na.rm=TRUE) # <-- ranks, not sensitivities!
                      , rankMeanOUT =   mean(rankNotTargetedVec, na.rm=TRUE) # <-- ranks, not sensitivities!
                      ,       medIN = median(sensTargetedVec   , na.rm=TRUE) # <-- growth (probably log10(exp/control): higher = more deadly drug)
                      ,      medOUT = median(sensNotTargetedVec, na.rm=TRUE) # <-- growth / sensitivity (higher = less growth / more deadly drug)
                      , sensMedAbsDev = theSensMedAbsDev
                      , rankMedAbsDev = theRankMedAbsDev
                      ,      meanIN =   mean(sensTargetedVec   , na.rm=TRUE) ## Mean of the sensitivity of this del-orf to the drugs that are annotated to target the pathway that we are currently looking at
                      ,     meanOUT = mean(sensNotTargetedVec  , na.rm=TRUE) ## Mean of the sensitivity of this del-orf to the drugs that are annotated to target some OTHER pathway than the current one
                      , overlap = overlapAmt
                      ,     orf = orf
                      ,    gene = theRealGeneName
                      ,   nTargetingModule  = length(sensTargetedVec) ## number of drugs targeting this module. Note that this is being inefficiently stored on a per-gene basis, instead of per-module
                      , nNotTargetingModule = length(sensNotTargetedVec)
                      ##, desc  <- gvYeastGeneDescriptions[orf, 3]
                      )
          
          ##t.test(sensTargeted, sensNotTargeted,alternative="greater")
     } # end "if"
     return(ttt)
} ## end function




agwSetUpDiffSensVsCentralityScoreGlobalVariables <- function() {
     ## set here are: 1. gvModVec
     ##               2. gvNeighborStatsList
     ##               2. gvNeighborStatsHash
     ##               4. gvMemberStatsList
     
     agwGlobalLoad("gvPathwayMembershipHash")
     agwGlobalLoad("gvGeneSLPathway")
     agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")
     agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
     agwGlobalLoad("gvDrugTargetMapping")
     agwGlobalLoad("gvYeastGeneNameHash")

     agwGlobalLoad("gv.HOM.RANKS.BY.DRUG") ## <-- used by agwFigureOutGeneInfo
     #unique(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,"TARGET_MODULE"])

     DELETION_GENE_INDEX = 1 ## index in gvGeneSLPathway
     MODULE_NAME_INDEX   = 2 ## index in gvGeneSLPathway

     if (!agwHasContent("gvModVec")) {
          gvModVec <<- list() ## <-- global in scope
          numRowsOfPredictions <- nrow(gvGeneSLPathway)
          for (i in 1:numRowsOfPredictions) {
               ##thisGene = gvGeneSLPathway[[DELETION_GENE_INDEX]]
               thisMod  = gvGeneSLPathway[i, MODULE_NAME_INDEX]
               if (is.null(thisMod) || (nchar(thisMod) == 0)) { print("Skipping empty module") ; next ; } ## skip this iteration of the "for" loop
               if (is.null(gvModVec[[thisMod]])) { gvModVec[[thisMod]] <<- vector() }
               gvModVec[[thisMod]] <<- append(gvModVec[[thisMod]], i) ## store the index of where to find this module's SL hits in the gvGeneSLPathway data structure
               if (i %% 25000 == 0) { print(paste("Populating gvModVec: Done with",i,"of", numRowsOfPredictions,"rows in the gvGeneSLPathway mutant-vs-pathway prediction list.")) }
          }
          print(paste("Status report: read in a total of this many pathways:", length(gvModVec)))
     }
     
     drugsThatWereTested = colnames(GLOBAL.HOM.SENS.COMPENDIUM)

     if (!agwHasContent("gvNeighborStatsList") || !agwHasContent("gvMemberStatsList") || !agwHasContent("gvNeighborStatsHash")) {
          gvNeighborStatsList <<- list() ## will save all the t-test / overlap info for any ORFs with an SL neighborhood with at least some overlap with the pathway in question
          gvMemberStatsList   <<- list() ## will save all the t-test / overlap info for ORFs that are actually IN this parthawy
          allOrfNamesTestedAnywhereVec <- rownames(GLOBAL.HOM.SENS.COMPENDIUM)
          agwGlobalLoad("gvDrugTargetMapping")
          gvNeighborStatsHash <<- agwNewHash(size=length(gvDrugTargetMapping))
          
          cat("Populating gvModVec with every ORF that is in the SL network AND is listed as the target of any drug (this takes QUITE A WHILE (~1 minute per drug), so be patient!)...")

          for (i in 1:length(gvDrugTargetMapping)) {
               ## This takes a very long time. If you want to speed it up, you can check only a few drugs if you want!
               #print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");print("DEBUGGING!!!!!!!!!!!!!!!!!!!!");
               
               mo <- names(gvDrugTargetMapping)[i] ## <-- a single module name
               cat("\n") ; cat(mo) # prints out the module name
               
               drugsThatTargetMo            <- names(gvDrugTargetMapping[[i]]) ## <-- a vector of *all* the drugs that are listed to target "mo"
               drugsThatDoNotTargetMo       <- setdiff(rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME), drugsThatTargetMo)
               testedDrugsThatTargetMo      <- intersect(drugsThatTargetMo, drugsThatWereTested)
               testedDrugsThatDoNotTargetMo <- intersect(drugsThatDoNotTargetMo, drugsThatWereTested)
               orfsInThisModuleList         <- agwHashGet(hash=gvPathwayMembershipHash, key=mo)
               orfInModuleHash              <- agwNewHash()
               agwHashPutMembership(orfInModuleHash, keys=orfsInThisModuleList, setValue=TRUE)
               
               ##agwPrint("Number of orfs in this module: ", length(orfInModuleHash))
               ##agwPrint(agwHashKeys(orfInModuleHash))
               ##browser()
               
               # ============================================================================================
               j <- 0
               for (memberORF in orfsInThisModuleList) {
                    singleORFInfo <- agwFigureOutGeneInfo(orf=memberORF, indexInSLCollection=NULL
                                                          , testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo)
                    if (!is.null(singleORFInfo)) {
                         if (j %% 10 == 1) { cat(":") } ## <-- "progress bar"
                         gvMemberStatsList[[mo]][[singleORFInfo$orf]] <<- singleORFInfo ## <-- note: these may be OVERWRITTEN later
                    }
                    j <- j+1
               }
               # ============================================================================================
               j <- 0
               for (geneIdx in gvModVec[[mo]]) {
                    ## Go through all the ORFs that are part of this module's SL neighborhood, and add their data to the neighbor stats.
                    ## Note that these are only the ORFs that actually *have* non-zero SL overlap.
                    ## The geneIdx runs through the various indices for ORFs that had significant overlap with the pathway in question
                    theOrfName    <- gvGeneSLPathway[geneIdx, DELETION_GENE_INDEX]
                    singleORFInfo <- agwFigureOutGeneInfo(orf=theOrfName, indexInSLCollection=geneIdx, testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo)
                    if (!is.null(singleORFInfo)) {
                         if (j %% 10 == 1) { cat(".") } ## <-- "progress bar"
                         gvNeighborStatsList[[mo]][[theOrfName]] <<- singleORFInfo
                         
                         agwHashOfHashesPut(gvNeighborStatsHash, key1=mo, key2=theOrfName, value=singleORFInfo)
                         
                         if (agwHashContains(orfInModuleHash, theOrfName)) {
                              gvMemberStatsList[[mo]][[theOrfName]] <<- singleORFInfo ## <-- overwrite it with the new info that actually contains the degree of SL overlap!
                              ## Line above: note that we are modifying the MEMBER stats list here! This is on purpose!
                         }
                    }
                    j <- (j + 1)
               } ## end of looping through all the modules with known target/indicators ORFs AND not too many targets
               # =======================================================
               j <- 0
               for (orf in allOrfNamesTestedAnywhereVec) {
                    ## Here, we handle the ORFs with zero-score SL overlap.
                    ## Note that we do not want to overwrite the non-zero SL overlap ORFs' entries, so if we find
                    ## that we already had data
                    
                    if (!is.null(agwHashOfHashesGet(hash=gvNeighborStatsHash, key1=mo, key2=orf))) {
                         next; ## Don't overwrite things here...
                    } else {
                         info <- agwFigureOutGeneInfo(orf=orf, indexInSLCollection=NULL,
                                                      testedDrugsThatTargetMo, testedDrugsThatDoNotTargetMo)
                         if (!is.null(info)) {
                              if (j %% 10 == 1) { cat("x") } ## <-- "progress bar"
                              gvNeighborStatsList[[mo]][[orf]] <<- info
                              agwHashOfHashesPut(gvNeighborStatsHash, key1=mo, key2=orf, value=singleORFInfo)
                         }
                    }
                    j <- (j+1)
               }
               # ============================================================================================
               ##agwPrint("Member: ", length(gvMemberStatsList)) ; agwPrint("Neighbor: ", length(gvNeighborStatsList))
          }
     }
     
     if (ADD.DEBUG.FAKE.DRUGS) { ## Adding a fake "good results, which should come to the top of the predictions list" drug for testing! The drug is named "FakeDrug," and there are 5 replicates.
          TESTER <- "TEST_MODULE"
          for (i in 1:20) {
               fakeORFName = agwGlue("FakeORF_", i)
               fakeGeneName = agwGlue("FAKE_", i)
               gvNeighborStatsList[[TESTER]][[fakeORFName]] <<- list(   meanIN = (i-3)*1.5
                                                                     , meanOUT = 4
                                                                     ,   medIN = (i-3)*1.5
                                                                     ,  medOUT = 4
                                                                     , overlap = i*0.7
                                                                     ,     orf = fakeORFName
                                                                     ,    gene = fakeGeneName
                                                                     , nTargetingModule    = 5
                                                                     , nNotTargetingModule = 5
                                                                     )
          }
          for (i in 1:5) {
               fakeDrugName <- agwGlue("fakeDrug_", i)
               gvDrugTargetMapping[[TESTER]][[fakeDrugName]] <<- 1
          }
          ##gvMemberStatsList[[TESTER]]
     }
     # If you want to see the distribution of any of the "neighborstatslist" stats:
     # WHICH_DRUG_IDX = 2  # <-- or whichever drug you want to check out
     #  hist(sapply(gvNeighborStatsList[[WHICH_DRUG_IDX]], function(a) { a$sensMedAbsDev }), xlab="Median Abs Deviation", col="brown")
     
} ## end of "agwSetUpDiffSensVsCentralityScoreGlobalVariables"


# =======================================================

# Calculating differential sensitivity:
## Suppose the 5 drugs that target "cell wall" super-duper kill this particular ORF:
## Then the orf growth values will be like:
##  4.8, 2.4, 4.8, 5.0, 9.1
## And the other drugs will look like:
##  0.4, 0.1, 0.1, 0.4, -0.1, -0.2, 0.4, 0.7, 1.2, -0.4  (negative means the del. ORF mutant grew *better* with the drug)
## So the means would be (say)
##  3.5 for the "drugs that target this pathway"
## and:
##  0.3 for the "drugs that do not target this pathway"
## Now if we subtract, we get:
## 3.5 - 0.3 = 3.2, which means that the ORFs that are in the pathway that is claimed to be targeted by these drugs
## are in fact more sensitive to those drugs than to other drugs

isDifferentiallySensitive <- function(sensAmt) {
     return (sensAmt > sensThresh)
     ## This is the current threshold for differential
     ## sensitivity i.e., "was the deletion ORF mutant in
     ## fact more sensitive to drugs that supposedly targeted
     ## the pathway in question, as compared to other drugs?"
}


# =======================================================

calcManyDiffSens <- function(dataList) { ## Calculate differential sensitivity (sens. in drugs that target this pathway vs sens for other drugs)
     # You call this on a LIST, where each element of the list
     # is yet ANOTHER list, of the type defined
     # The data for "dataList" is set in agwDrugCalcFunctions. Search for "ttt" or for "medIN" .
     calcSingleDiffSens <- function(aList) {
     # This is differential sensitivity as defined in the paper.
     # Note that we divide the median difference by the
     # median absolute deviation ("mad(...)"), as calculated
          if (!is.null(aList) && !is.na(aList)
              && (length(aList$medIN) > 0)  && !is.null(aList$medIn)  && !is.na(aList$medIn)
              && (length(aList$medOUT) > 0) && !is.null(aList$medOut) && !is.na(aList$medOut)
              && (length(aList$sensMedAbsDev) > 0)
              && (aList$sensMedAbsDev > 0)) {
               ## This is how diff. sens is defined in the paper: DIS(g, k)
               ## (median sens of hitting drugs - median sens of other drugs)
               ## -------------divided by--------------------------
               ## (median absolute deviation of all sens. for this ORF)
               return((aList$medIN - aList$medOUT)/(aList$sensMedAbsDev))  ## <-- differential sensitivity is defined here!
               # return(a$medIN - a$medOUT) # <-- old-style
          } else {
               return(NA) # No valid differential sens. Maybe there were no drugs in the "hitting this pathway" group.
          }
     }
     
     z <- sapply(dataList, function(a) { calcSingleDiffSens(a) })
     if (length(z) > 0) { return(z) }
     else { return(c(NA)) }
}



## <> MAIN <>
## ========== MAIN ===============


