# source("/Users/alexgw/R/src/DrugEval.R"); # <-- to reload this file


source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");


source("/Users/alexgw/R/src/agwTEST.R"); # <-- all the testing stuff!!

#agwHandleDifferentialSensitivity()



testpred <- agwReadFileIntoDataFrame(agwInputPredictionFilePath(INPUT.DATA.ROOT.DIR, "GoMyers", "Ptwy", "HOM", "neigh", "wilcoxon")) ; # <-- the input file's complete path

# Suppose you want a histogram of the sensitivities of one *particular* drug in the (say) HOM compendium:
# hist(ztab.hom.compendium[["Hillenmeyer/papuamide-B_0.7_ug#1"]],col="lightblue",breaks=60); agwLine(0, 'v')

# In order to get the entire list of drugs, try colnames(ztab.hom.compendium)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&& START HANDLING PLOTS &&&&&&&&&&&&&&&&&
# &&&&&&
# &&&&
# &&&
# &&
# &
agwHandlePlots <- function(modset, modtype, zygosity, distance, scoretype) {
	#
	cat("\n\n***\n**\n*\n\nDEBUGGING & TESTING ENABLED. THIS OVERWRITES THE INPUT ARGUMENTS TO agwHandlePlots. COMMENT THIS OUT IN THE REAL RUN!!\n\n*\n**\n***\n\n") ; modset = "GoMyers" ; modtype = "Ptwy" ; zygosity = "HOM"; distance = "neigh"; scoretype="wilcoxon" ; # <-- for testing only! comment this out later
	
	infile <- agwInputPredictionFilePath(INPUT.DATA.ROOT.DIR, modset, modtype, zygosity, distance, scoretype) ; # <-- the input file's complete path
	print(agwGlue("Reading from file <",infile,">"));
		
	#SOURCE.PRED.FILE  = "~/R/Data/Pred/GoMyers/Ptwy/HOM/neigh/Pred/wilcoxon_enrichment_pval.tab" # PRED.FILES[1] #
	
	pred <- agwReadFileIntoMatrix(infile)
	
	# ================
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
	
	agwPredBoxPlot(tplot.mat.sorted.by.max.transposed
		, agwOutPath(modset, modtype, zygosity, distance, scoretype, "Predictions_by_module_SortedByMax.png")
		, numModulesTested=ncol(tplot.mat.sorted.by.max.transposed))
	
	agwPredBoxPlot(tplot.mat.alpha.transposed
		, agwOutPath(modset, modtype, zygosity, distance, scoretype, "Predictions_by_module_SortedAlphabetically.png")
		, numModulesTested=ncol(tplot.mat.alpha.transposed))
	
	agwPredBoxPlot(tplot.mat.sorted.by.max
		, agwOutPath(modset, modtype, zygosity, distance, scoretype, "Predictions_by_drug_SortedByMax.png")
		)
	agwPredBoxPlot(tplot.mat.alphabetical
		, agwOutPath(modset, modtype, zygosity, distance, scoretype, "Predictions_by_drug_SortedAlphabetically.png")
		)

}
# &
# &&
# &&&
# &&&&
# &&&&&&
# &&&&&&&&&&&&&&&&&& DONE HANDLING PLOTS &&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


agwHandlePlots("GoMyers","Ptwy","HOM","neigh","wilcoxon") # <-- for debugging!

agwHandlePlots("GoMyers","Ptwy","HOM","direct","wilcoxon") # <-- for debugging!
agwHandlePlots("GoMyers","Ptwy","HET","direct","wilcoxon") # <-- for debugging!


PRED.MODULES = c("Compendium","GoMyers")
PRED.MODTYPE = c("Gene","Ptwy")
PRED.ZYGOSITY = c("HET","HOM")
PRED.DISTANCE = c("direct","neigh")
PRED.SCORETYPE = c("wilcoxon","ttest","ks","gsea")

globalMod   = NULL
globalType  = NULL
globalZy    = NULL
globalDist  = NULL
globalScore = NULL

for (globalMod in PRED.MODULES) {
	for (globalType in PRED.MODTYPE) {
		for (globalZy in PRED.ZYGOSITY) {
			for (globalDist in PRED.DISTANCE) {
				for (globalScore in PRED.SCORETYPE) {
					
					#agwHandlePlots(globalMod, globalType, globalZy, globalDist, globalScore)
					agwHandlePredictionEvalPlot(globalMod, globalType, globalZy, globalDist, globalScore)
					
				}	
			}
			
		}		
	}
}

# Summarize the various properties of the positive control pathways.
# Note that this only needs to be done once--it is NOT on a per-scoretype basis
agwHandlePositiveControlGlobalSummary()

agwPlotSensitivitySummaries("HOM")
agwPlotSensitivitySummaries("HET")

# =================================================================

agwHandlePredictionEvalPlot("GoMyers","Ptwy","HOM","neigh","wilcoxon")

# =================================================================

#agwGlue(GLOBAL.OUTPUT.PATH, ZYGOSITY, "_Hisogram_of_Mean_Sensitivities",".png")



#stripchart(alz[1:10], method = 'jitter',pch='.')


#par(srt=90)
#text(x=seq(1, 4, by=1), y=rep(-5,4), labels=colnames(plotz), pos =1, xpd = TRUE) 
#mtext(colnames(boxplotz), at=seq(1,ncol(plotz)), side=2, outer=FALSE, adj=1,padj=1,las=LAS.REGULAR.TEXT, line=1)



#den <- density(log(positives))
#den$x <- exp(den$x)
#plot(den, log="x", main="something meaningful")
#c(0,1,2,3,4,5,6,7,8,9,10,100,max(positives[,5]))


