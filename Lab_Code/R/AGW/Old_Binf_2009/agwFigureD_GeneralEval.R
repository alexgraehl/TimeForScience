

#    source("/Users/alexgw/R/src/agwFigureD_GeneralEval.R");    # <-- to reload this file

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

#source("/Users/alexgw/R/src/agwTEST.R"); # <-- all the testing stuff!!

#agwHandleDifferentialSensitivity()

#testpred <- agwReadFileIntoDataFrame(agwInputPredictionFilePath(DATA_ROOT_DIR, "GoMyers", "Ptwy", "HOM", "neigh", "wilcoxon")) ; # <-- the input file's complete path

# Suppose you want a histogram of the sensitivities of one *particular* drug in the (say) HOM compendium:
# hist(ztab.hom.compendium[["Hillenmeyer/papuamide-B_0.7_ug#1"]],col="lightblue",breaks=60); agwLine(0, 'v')

# In order to get the entire list of drugs, try colnames(ztab.hom.compendium)

FIGURE_D_SUBDIR     = "FigureD_General"
FIGURE_D2_SUBDIR    = "FigureD2_PredSummary"

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&& START HANDLING PLOTS &&&&&&&&&&&&&&&&&
# &&&&&&
# &&&&
# &&&
# &&
# &

## =================================================================

agwInputPredictionFilePath <- function(root, modset, modtype, zygosity, distance, scoretype) {
     UNIX.DIRECTORY.CHAR = '/';
     return(paste(root, PREDICTIONS_DIRNAME_IN_ROOT, modset, modtype, zygosity, distance, "Pred", agwGlue(scoretype, "_enrichment_pval.tab"), sep=UNIX.DIRECTORY.CHAR) );
}



## =================================================================
agwHandlePositiveControlGlobalSummary <- function() {
     ## Here, let's analyze the positive controls.
     ## First, we load in the entire table of positive controls (known targets)

     agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")
     agwGlobalLoad("GLOBAL.POSITIVES.NO.DUPS.DATAFRAME")

     ## GLOBAL.POSITIVES.NO.DUPS.DATAFRAME is a similar table, but without separate entries 
     ## for each drug dosage. It's much shorter, and the columns are reordered
     COLUMN.WITH.NUM.TARGET.ORFS = 5

     max.number.of.target.orfs.for.the.positive.control.pathway <- 300 ## <-- this is somewhat arbitrary

     posWithoutTooManyORFTargets <- which(GLOBAL.POSITIVES.NO.DUPS.DATAFRAME[,COLUMN.WITH.NUM.TARGET.ORFS] <= max.number.of.target.orfs.for.the.positive.control.pathway)

     ## OUTPUT PLOT: A two-part histogram showing the
     ## distribution of the count of annotated ORF targets for each drug with a known target.

     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D2_SUBDIR)
                    , file="PositiveControls_Number_of_Target_ORFs"
                    , pointsize=HIST.DEFAULT.POINTSIZE*0.75
                    , width=(hist.default.height*2)
                    , height=hist.default.height)
     split.screen(c(1,2))
     screen(1) ; agwHistogram(posWithoutTooManyORFTargets, out=NULL, breaks=10, xlab="No. target ORFs for a drug (>=50 only)")
     screen(2) ; agwHistogram(as.numeric(GLOBAL.POSITIVES.NO.DUPS.DATAFRAME[,COLUMN.WITH.NUM.TARGET.ORFS]), 	out=NULL, breaks=50, xlab="No. target ORFs for a drug")
     close.screen(all=TRUE)
     agwFinishPlot()

}
## =================================================================


## =================================================================
agwHandlePredictionEvalPlot <- function(modset, modtype, zygosity, distance, scoretype) {

     infile <- agwInputPredictionFilePath(DATA_ROOT_DIR, modset, modtype, zygosity, distance, scoretype)
     agwPrint("agwHandlePredictionEvalPlot: reading from the input file ", infile)

     ##thepred <- agwReadFileIntoDataFrame("~/R/Data/Pred/GoMyers/Ptwy/HOM/neigh/Pred/wilcoxon_enrichment_pval.tab")

     thepred <- agwReadFileIntoDataFrame(infile)

     agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")

     ## "mi" is a matrix of data about which positive control target pathway is associated with which drug
     ## It also contains a bit of extra data, like how many ORF targets are in the supposed target set
     ## Finds the drug names AND target pathways for each name
     mi <- cbind(rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME)
                 , agwGlue("SL_Neighbors_of:",GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,"TARGET_MODULE"])
                 , GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,"COUNT"])

     ## note that above, we add the string "SL_Neighbors_of:" to the known-target indication. That's because the genes that we expect to be enriched for being-killed-ness are the ones that are the SL neighbors of the "true" target.
     
     stopifnot(length(rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME)) == length(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,"TARGET_MODULE"]))
     stopifnot(length(rownames(GLOBAL.POSITIVES.COMPLETE.DATAFRAME)) == length(GLOBAL.POSITIVES.COMPLETE.DATAFRAME[,"COUNT"]))

     ##mi <- mi[1:40,] ## subsetting for debugging/testing...

     pos.performance <- rep(NA,nrow(mi))   ## How well the "positive control" known targets did in being predicted. Initialize them to "na" to start

     DRUGNAME_COL_IN_MI_IDX         = 1
     DRUG_TARGET_PTWY_COL_IN_MI_IDX = 2 ## <-- a pathway

     for (i in seq(1,nrow(mi)) ) {
          drug   <- mi[i, DRUGNAME_COL_IN_MI_IDX]
          target <- mi[i, DRUG_TARGET_PTWY_COL_IN_MI_IDX]
          val    <- thepred[target, drug]
          if (!is.null(val)) { pos.performance[i] <- val }
     }

     if (all(is.na(pos.performance))) {
          print("WARNING: All of pos.performance was NA for some reason... probably you are trying to run this on GENE data. But all the data is for pathways. Sorry.")
          return()
          #browser()
     }
     
     numMultipleTests <- nrow(thepred)
     confidence.0.05  <- (-log10(0.05 / numMultipleTests))
     confidence.0.01  <- (-log10(0.01 / numMultipleTests))
     confidence.0.001 <- (-log10(0.001 /numMultipleTests))

     pos.performance[pos.performance < ENRICH.LIMIT.SMALL.P.VALUE] <- ENRICH.LIMIT.SMALL.P.VALUE ## anything LESS than the minimum value gets clipped (avoids erroneous p-values of 0 due to rounding)
     ##boxplot.row.height*ncol(pos.performance)
     
     longLabels     <- agwGlue(mi[,1], " --> ",mi[,2]," (c = ",mi[,3],")")
     abridgedLabels <- agwAbridgeModuleNames(longLabels)
     
     agwPrint("Showing the performance of the positive control pathway performance...")

     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D2_SUBDIR)
                    , file=paste(modset, modtype, zygosity, distance, scoretype, "TargetPathwayPerformanceOnly", sep="_")
                    , pointsize=DEFAULT.POINTSIZE
                    , width=2*STANDARD.PLOT.WIDTH
                    , height=(DEFAULT.ROW.HEIGHT * length(pos.performance)))
     par(mar=c(boxplot.m.bottom, boxplot.m.left, boxplot.m.top, boxplot.m.right))
     
     pos.performance.colors = ifelse(is.na(pos.performance), "dark gray", "black")
     
     dotchart(-log10(pos.performance), labels=abridgedLabels, pch=PCH.DIAMOND, bg="red", cex=1.0, lcolor="gray", xlab="-log10(Pval) of drug's annotated target", col=pos.performance.colors)
     agwLine(confidence.0.05,'v') ; agwLine(confidence.0.01,'v') ; agwLine(confidence.0.001,'v')
     agwFinishPlot()


     print("Showing the positive control pathway scores vs the size of those pathways...")
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D2_SUBDIR)
                    , file=paste(modset, modtype, zygosity, distance, scoretype, "Positives--score_vs_size", sep="_")
                    , res=140
                    , width=1400
                    , height=1400)
     m.bottom = 5 ;  m.left = 5 ; m.top = 4 ; m.right = 4;
     par(mar=c(m.bottom, m.left, m.top, m.right))
     plot(mi[,3],-log10(pos.performance),pch=PCH.DIAMOND,bg="green", xlab="Module Size (No. ORFs)", ylab="Target predictions for this module (-log10(p): higher = better).", main="Positive Controls Only: Module scores & module sizes")
     agwLine(confidence.0.05,'h') ; agwLine(confidence.0.01,'h') ; agwLine(confidence.0.001,'h')
     text(confidence.0.05 + 2.15, pos=POS.RIGHT,  "Bonferroni-corrected p=0.05")
     text(confidence.0.01 + 0.30, pos=POS.RIGHT,  "Bonferroni-corrected p=0.01")
     text(confidence.0.001 - 2.3, pos=POS.RIGHT, "Bonferroni-corrected p=0.001")
     agwFinishPlot()

     ## Double checking (spot check): 0.0239 is the prediction value in "~/R/Data/Pred/GoMyers/Ptwy/HOM/neigh/Pred/wilcoxon_enrichment_pval.tab" for the pair "Hillenmeyer/floxuridine_3.37_um##2" and its target, "OUR_KNOWN_TARGET:DRUGBANK_MAPPED_FROM_HUMAN:floxuridine"
     ##thepred["OUR_KNOWN_TARGET:HAND_CURATED:protein-tyrosine_kinase","Hillenmeyer/pH8##6"]		
}

## =================================================================
## =================================================================
## =================================================================


agwHandlePredictionPlotsByDrug <- function(modset, modtype, zygosity, distance, scoretype) {
     #
     #cat("\n\n***\n**\n*\n\nDEBUGGING & TESTING ENABLED. THIS OVERWRITES THE INPUT ARGUMENTS TO agwHandlePredictionPlotsByDrug. COMMENT THIS OUT IN THE REAL RUN!!\n\n*\n**\n***\n\n") ; modset = "GoMyers" ; modtype = "Ptwy" ; zygosity = "HOM"; distance = "neigh"; scoretype="wilcoxon" ; # <-- for testing only! comment this out later
     
     infile <- agwInputPredictionFilePath(DATA_ROOT_DIR, modset, modtype, zygosity, distance, scoretype) ; # <-- the input file's complete path
     agwPrint("Reading from the file of predictions <",infile,">")
     #SOURCE.PRED.FILE  = "~/R/Data/Pred/GoMyers/Ptwy/HOM/neigh/Pred/wilcoxon_enrichment_pval.tab" # PRED.FILES[1] #
     pred <- agwReadFileIntoMatrix(infile)
     
     # Make two giant box plots showing the various predictions that were made for the various drugs
     # NAs are set to the worst P-value. P-values that are too close to zero (such that enrich.R reports them as zero) are set to the "ENRICH.LIMIT.SMALL.P.VALUE"
     tplot.mat <- pred
     tplot.mat[is.na(tplot.mat)] = WORST.P.VALUE # <-- all the NA values are set to the worst p-value (i.e., 1.0)
     tplot.mat[tplot.mat < ENRICH.LIMIT.SMALL.P.VALUE] = ENRICH.LIMIT.SMALL.P.VALUE # <-- anything that is under the "this p-value looks like zero" limit is set to some small p-value. P-values can't be zero, but sometimes enrich.R sends us 0 p-values anyway
     tplot.mat <- (-log10(tplot.mat)) # <-- now we log-scale the whole thing to make it easier to understand. Log base 10.
     
     #tplot.mat[tplot.mat < 5] = NA # <-- just for looking at only the extremely-confident values
     
     tplot.mat.sorted.by.max <- tplot.mat[,order(apply(tplot.mat, BY.COL, max, na.rm=TRUE))] # Order by the maximum values
     tplot.mat.alphabetical  <- tplot.mat[,order(colnames(tplot.mat), decreasing=TRUE)] # Order alphabetically by column name
     
     tplot.mat.sorted.by.max.transposed <- t(tplot.mat)
     tplot.mat.sorted.by.max.transposed <- tplot.mat.sorted.by.max.transposed[,order( apply(tplot.mat.sorted.by.max.transposed, BY.COL, max, na.rm=TRUE ))]
     
     tplot.mat.alpha.transposed <- t(tplot.mat.alphabetical)
     tplot.mat.alpha.transposed <- tplot.mat.alpha.transposed[,order(colnames(tplot.mat.alpha.transposed), decreasing=TRUE)]

     agwPrint("About to generate boxplot #1 of 4 for  <",infile,">... (this is slow!)")     
     agwPredBoxPlot(tplot.mat.sorted.by.max.transposed
                    , outFilename=agwGlobalOutPath(FIGURE_D2_SUBDIR, "/", modset, "_", modtype, "_", zygosity, "_", distance, "_", scoretype, "_Predictions_by_module_SortedByMax.png")
                    , numModulesTested=ncol(tplot.mat.sorted.by.max.transposed))

     agwPrint("About to generate boxplot #2 of 4 for  <",infile,">... (this is slow!)")
     agwPredBoxPlot(tplot.mat.alpha.transposed
                    , outFilename=agwGlobalOutPath(FIGURE_D2_SUBDIR, "/", modset, "_", modtype, "_", zygosity, "_", distance, "_", scoretype, "_Predictions_by_module_SortedAlphabetically.png")
                    , numModulesTested=ncol(tplot.mat.alpha.transposed))

     agwPrint("About to generate boxplot #3 of 4 for  <",infile,">... (this is slow!)")
     agwPredBoxPlot(tplot.mat.sorted.by.max
                    , outFilename=agwGlobalOutPath(FIGURE_D2_SUBDIR, "/", modset, "_", modtype, "_", zygosity, "_", distance, "_", scoretype, "_Predictions_by_drug_SortedByMax.png"))

     agwPrint("About to generate boxplot #4 of 4 for  <",infile,">... (this is slow!)")
     agwPredBoxPlot(tplot.mat.alphabetical
                    , outFilename=agwGlobalOutPath(FIGURE_D2_SUBDIR, "/", modset, "_", modtype, "_", zygosity, "_", distance, "_", scoretype, "_Predictions_by_drug_SortedAlphabetically.png"))

}


## =================================================================
## =================================================================
## =================================================================

agwPlotSensitivitySummaries <- function(ZYGOSITY=NULL) {

     agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
     agwGlobalLoad("GLOBAL.HET.SENS.COMPENDIUM")

     ## Suppose you want a histogram of the sensitivities of one *particular* drug in the (say) HOM compendium:
     ## hist(GLOBAL.HOM.SENS.COMPENDIUM[["Hillenmeyer/papuamide-B_0.7_ug##1"]],col="lightblue",breaks=60); agwLine(0, 'v')

     ## In order to get the entire list of drugs, try colnames(GLOBAL.HOM.SENS.COMPENDIUM)

     SENS.COLOR   = "pink" ## <-- unset! Should indicate a warning / error
     MEAN.COLOR   = "pink"
     MEDIAN.COLOR = "pink"

     if (ZYGOSITY == "HET") {
          SENS.COLOR = "orange"
          MEAN.COLOR = "red"
          MEDIAN.COLOR = "purple"
          WHICH.DATA.TABLE = GLOBAL.HET.SENS.COMPENDIUM
          
     } else if (ZYGOSITY == "HOM") {
          SENS.COLOR = "blue"	
          MEAN.COLOR = "dark green"
          MEDIAN.COLOR = "dark blue"
          WHICH.DATA.TABLE = GLOBAL.HOM.SENS.COMPENDIUM
     } else {
          stop("ERROR: ZYGOSITY was an unknown value. Error 123456654321.")
     }

     MEAN.SYMBOL   = 20 ## <-- pch item #20, not actually the string "20"
     MEDIAN.SYMBOL = '+'
     SENSITIVITY.TYPE  = "log_10(ctrl_growth/experiment_growth)"
     SENSITIVITY.LABEL = paste(SENSITIVITY.TYPE, "(higher = less growth when exposed to drug)")
     MEAN.SENSITIVITY.LABEL   = paste("Mean of", SENSITIVITY.LABEL)
     MEDIAN.SENSITIVITY.LABEL = paste("Median of", SENSITIVITY.LABEL)
     MAIN.LABEL        = agwGlue("mutant drug sensitivity (", ZYGOSITY, ")")
     MEAN.MAIN.LABEL   = paste("Mean",   MAIN.LABEL)
     MEDIAN.MAIN.LABEL = paste("Median", MAIN.LABEL)

     alz <- rbind(WHICH.DATA.TABLE)[,order(colnames(WHICH.DATA.TABLE))] ## Alphabetial sort by drug name (columns)

     ## Calculate the means and medians
     ## ================================
     hmean <- mean(alz, na.rm=TRUE) ## <-- Calculate MEAN values
     hmedian = apply(alz, BY.COL, median, na.rm=TRUE) ## <-- Calculate MEDIAN values

     mediansByMeanSort <- rbind(hmedian)[order(hmean)]
     medianNamesByMeanSort <- rbind(names(hmedian))[order(hmean)]
     ## ================================

     ## ================================
     ## Plot the MEANS and MEDIANS of the drug sensitivities, sorted by mean value

     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Mean_And_Median_Drug_Sensitivities_Sorted_by_Mean_Value"), pointsize=12,width=STANDARD.PLOT.WIDTH,height=STANDARD.PLOT.HEIGHT)
     plot(sort(hmean),ylab=MEAN.SENSITIVITY.LABEL,main=paste("Mean (O) and Median (+)", MAIN.LABEL)
          , xlab="Drug, mean value (higher=more sensitive)"
          , col=MEAN.COLOR
          , pch=MEAN.SYMBOL
          , ylim=c(min(hmean,hmedian,na.rm=TRUE),max(hmean,hmedian,na.rm=TRUE)))
     agwLine(0, 'h')
     ##points(sort(hmedian))
     points(mediansByMeanSort, col=MEDIAN.COLOR,pch=MEDIAN.SYMBOL)
     legend("bottomright", c("Mean","Median"), pch=c("O",MEDIAN.SYMBOL), title="Key",inset=0.025,col=c(MEAN.COLOR,MEDIAN.COLOR) )
     ##text(sort(hmean),names(sort(hmean)),pos=POS.RIGHT) ## Plot the names for the MEAN points
     ##text(mediansByMeanSort, medianNamesByMeanSort,pos=POS.RIGHT)
     ## The super-low guy is 7C18_50_uM (display the super high/low data points with: sort(hmean))
     agwFinishPlot()
     ## ================

     ## ================
     ## Plot the MEDIANS of all drug sensitivities, sorted by value
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Median_Drug_Sens_by_value"))
     plot(sort(hmedian),ylab=MEDIAN.SENSITIVITY.LABEL,main=MEDIAN.MAIN.LABEL
          , xlab="Drugs sorted by median value (higher=more sensitive)"
          , col=MEDIAN.COLOR
          , pch=MEDIAN.SYMBOL)
     agwLine(0, 'h')
     agwFinishPlot()
     ## ================

     ## ================
     ## Plot the MEDIANS and MEANS of the drug sensitivities, sorted ALPHABETICALLY

     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Mean_And_Median_Drug_Sensitivities_Sorted_Alphabetically")
                    , pointsize=12
                    , width=STANDARD.PLOT.WIDTH
                    , height=STANDARD.PLOT.HEIGHT)
     plot(hmean
          , main=MEAN.MAIN.LABEL
          , ylab=MEAN.SENSITIVITY.LABEL
          , xlab="Index of drug, sorted alphabetically"
          , col=MEAN.COLOR
          , pch=MEAN.SYMBOL)
     points(hmedian, col=MEDIAN.COLOR, pch=MEDIAN.SYMBOL)
     legend("bottomright", c("Mean","Median"), pch=c("O",MEDIAN.SYMBOL), title="Key", inset=0.025, col=c(MEAN.COLOR,MEDIAN.COLOR) )
     agwLine(0, 'h')
     agwFinishPlot()
     ## ================
     sort(hmean)
     ## ----------------


     ## ================
     ## Plot just the MEDIANS of all drug sensitivities, sorted ALPHABETICALLY
     ##FILE = agwGlue(GLOBAL.OUTPUT.PATH, ZYGOSITY, "_Medians_of_drug_sens_alphabetical",".png")
     ##print(FILE) ; png(FILE,pointsize=12,width=STANDARD.PLOT.WIDTH,height=STANDARD.PLOT.HEIGHT)
     #plot(hmedian,ylab=MEDIAN.SENSITIVITY.LABEL,main=MEDIAN.MAIN.LABEL,xlab="Index of drug, sorted  alphabetically",col=MEDIAN.COLOR, pch=MEDIAN.SYMBOL)
     #agwLine(0, 'h')
     ##legend("bottomright", c("Mean","Median"), pch=c("O",MEDIAN.SYMBOL), title="Key",inset=0.025,col=c(MEAN.COLOR,MEDIAN.COLOR))

     ## ================
     sort(hmedian)
     ## ----------------

     ## ================
     ## Calculate how many items in alz were NOT NA. i.e., how many mutants were tested for each drug

     
     alz.notNA.count <- apply(alz, BY.COL, function(a) { length(na.omit(a))})
     
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Histogram_of_deletion_mutant_stats"), pointsize=12,width=STANDARD.PLOT.WIDTH,height=STANDARD.PLOT.HEIGHT)
     hist(alz.notNA.count,breaks=80, col=SENS.COLOR,labels=TRUE,xlab="Number of deletion mutants tested.\nYeast has ~5000 nonessential genes and ~1000 essential genes.",ylab= agwGlue("Count of drug experiments (Total: ",length(alz.notNA.count)," drugs tested)"),main=paste(ZYGOSITY,": The number of gene deletion mutants tested in\nvarious drug+mutant experiments."))
     agwFinishPlot()
     
     ## To show a list of the number of deletion mutants tested in each study:
     sort(alz.notNA.count)

     ## The 14 items in the Rcy1 group have the most (5927) genes tested. This presumably includes both essential and nonessential genes. Rcy1/C700-1459_5.625uM_08_11_27_Lokey_t06
     ## Everything else is around 4000-4750, which is more plausible for only the non-essential genes.
     ## ================

     sens.all.hist.title = agwGlue("sensitivities for ",ncol(alz)," tested drugs (log10(ctrl/test), higher=less growth)")
     sens.all.hist.xlab  = agwGlue(" growth inhibition, log10 scale.\n1.0 = 10-fold less growth with drug, vs. control.")
     ## ================
     ## Histogram of the MEAN sensitivities
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR)
                    , file=agwGlue(ZYGOSITY, "_histogram_of_mean_sensitivities")
                    , pointsize=12, width=STANDARD.PLOT.WIDTH, height=STANDARD.PLOT.HEIGHT)
     hist(hmean, breaks=30, col=MEAN.COLOR, labels=TRUE, xlab=agwGlue("Mean ", sens.all.hist.xlab)
          , main=paste(ZYGOSITY, "-- Mean", sens.all.hist.title))
     rug(hmean,col=MEAN.COLOR)
     agwLine(0, 'v')
     agwFinishPlot()
     ## ================

     ## ================
     ## Histogram of the MEDIAN sensitivities
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR)
                    , file=agwGlue(ZYGOSITY, "_Histogram_of_median_sensitivities")
                    , pointsize=12, width=STANDARD.PLOT.WIDTH, height=STANDARD.PLOT.HEIGHT)
     hist(hmedian, breaks=20, col=MEDIAN.COLOR, labels=TRUE, xlab=agwGlue("Median ", sens.all.hist.xlab)
          , main=paste(ZYGOSITY, "-- Median", sens.all.hist.title))
     rug(hmedian,col=MEDIAN.COLOR)
     agwLine(0, 'v')
     agwFinishPlot()
     ## ================


     ## ================
     ## Output to a file with a plot of the mean sensitivities of each drug, in alphabetical order.
     ## Plots a dot chart of drug sensitivity.
     ## the vector to plot in dotchart is expected to be a vector of means or medians, with labels for each drug.

     agwDotCharter <- function(plotItVec, main, bg, xlab, type) {
          agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR)
                         , file=agwGlue(ZYGOSITY, "_DotChart_Sensitivities_", type, "_Alphabetical")
                         , width=1200
                         , height=(100 + 20*length(plotItVec)) );
          dotchart(plotItVec, main=main, pch=PCH.DIAMOND, bg=bg, cex=1.0, lcolor="gray", xlab=xlab)
          abline(v=seq(from=-0.5,to=0.5,by=0.25), lty=1, lwd=1,col="orange")
          agwLine(0, 'v')
          agwFinishPlot()
     }

     agwDotCharter(hmean  , main=MEAN.MAIN.LABEL  , bg=MEAN.COLOR  , xlab=MEAN.SENSITIVITY.LABEL  , type="MEAN")
     agwDotCharter(hmedian, main=MEDIAN.MAIN.LABEL, bg=MEDIAN.COLOR, xlab=MEDIAN.SENSITIVITY.LABEL, type="MEDIAN")
     ## ================

      ## Variables used for the box plot
     boxplotz <- alz ##[,1:20]

     ##pos.performance.colors = pos.performance
     ##pos.performance.colors[!is.na(pos.performance)] = "black"
     ##pos.performance.colors[is.na(pos.performance)] = "dark gray"
     ## Make a giant box plot showing the sensitivities of everything
     ## You have to run the variables above first, since they get used here
     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Boxplot_sensitivities_summary_alphabetical"), pointsize=12,width=boxplot.width, height=100+boxplot.row.height*length(boxplotz))
     par(mar=c(boxplot.m.bottom, boxplot.m.left, boxplot.m.top, boxplot.m.right))
     boxplot(boxplotz, horizontal=TRUE, xlab=SENSITIVITY.LABEL, main=paste(ZYGOSITY, "box plot of drug sensitivity.\nEach point is a single drug+mutant experiment."), pch=21,bg=SENS.COLOR,col=SENS.COLOR,axes=TRUE,las=LAS.HORIZONTAL.TEXT)
     agwLine(0, 'v')
     abline(v=c(-1,1), lty=LTY.DASHED.LINE,col="orange", lwd=2)
     agwFinishPlot()
     ## ================

     ## ================
     ## Make a giant box plot showing the sensitivities of everything
     ## You have to run the variables above first, since they get used here
     boxplotzValueSorted <- rbind(boxplotz)[,order(hmean)]

     agwPreparePlot(directory=agwGlobalOutPath(FIGURE_D_SUBDIR), file=agwGlue(ZYGOSITY, "_Boxplot_sensitivities_summary_sorted_by_mean_value"), pointsize=12,width=boxplot.width, height=100+boxplot.row.height*length(boxplotz))
     par(mar=c(boxplot.m.bottom, boxplot.m.left, boxplot.m.top, boxplot.m.right))
     boxplot(boxplotzValueSorted, horizontal=TRUE, xlab=SENSITIVITY.LABEL, main=paste(ZYGOSITY, "box plot of drug sensitivity.\nEach point is a single drug+mutant experiment."), pch=21,bg=SENS.COLOR,col=SENS.COLOR,axes=TRUE,las=LAS.HORIZONTAL.TEXT)
     agwLine(0, 'v')
     abline(v=c(-1,1), lty=LTY.DASHED.LINE,col="orange", lwd=2)
     agwFinishPlot()
     ## ================

}


## =================================================================


# &
# &&
# &&&
# &&&&
# &&&&&&
# &&&&&&&&&&&&&&&&&& DONE HANDLING PLOTS &&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


## =================================================================


agwPlotPredEvalPlots <- function() {
     PRED.MODULES = c("Compendium","GoMyers")
     PRED.MODTYPE = c("Gene","Ptwy")
     PRED.ZYGOSITY = c("HET","HOM")
     PRED.DISTANCE = c("direct","neigh")
     PRED.SCORETYPE = c("wilcoxon") #,"ttest","ks","gsea")

     ##globalMod   = NULL ; globalType  = NULL ; globalZy    = NULL ; globalDist  = NULL ; globalScore = NULL

     for (aMod in PRED.MODULES) {
          for (aType in PRED.MODTYPE) {
               for (aZy in PRED.ZYGOSITY) {
                    for (aDist in PRED.DISTANCE) {
                         if (aZy == "HET" && aDist == "neigh") {
                              next; # skip these for now... I think this is actually technically included in the direct results
                         }
                         
                         for (aScore in PRED.SCORETYPE) {
                              agwHandlePredictionPlotsByDrug(aMod, aType, aZy, aDist, aScore)
                              agwHandlePredictionEvalPlot(aMod, aType, aZy, aDist, aScore)
                              
                         }	
                    }
                    
               }		
          }
     }
} ## <-- end function

# Summarize the various properties of the positive control pathways.
# Note that this only needs to be done once--it is NOT on a per-scoretype basis

agwPlotSensitivitySummaries("HOM") ## goes into FIGURE_D_SUBDIR
agwPlotSensitivitySummaries("HET") ## goes into FIGURE_D_SUBDIR
agwHandlePositiveControlGlobalSummary() ## goes into FIGURE_D2_SUBDIR
agwPlotPredEvalPlots() ## Goes into FIGURE_D2_SUBDIR

#agwHandlePredictionPlotsByDrug("GoMyers", "Ptwy", "HOM", "neigh", "wilcoxon")

# =================================================================

#agwHandlePredictionEvalPlot("GoMyers","Ptwy","HOM","neigh","wilcoxon")

# =================================================================

#agwGlue(GLOBAL.OUTPUT.PATH, ZYGOSITY, "_Histogram_of_Mean_Sensitivities",".png")

#stripchart(alz[1:10], method = 'jitter',pch='.')

#par(srt=90)
#text(x=seq(1, 4, by=1), y=rep(-5,4), labels=colnames(plotz), pos =1, xpd = TRUE) 
#mtext(colnames(boxplotz), at=seq(1,ncol(plotz)), side=2, outer=FALSE, adj=1,padj=1,las=LAS.REGULAR.TEXT, line=1)

#den <- density(log(positives))
#den$x <- exp(den$x)
#plot(den, log="x", main="something meaningful")
#c(0,1,2,3,4,5,6,7,8,9,10,100,max(positives[,5]))





