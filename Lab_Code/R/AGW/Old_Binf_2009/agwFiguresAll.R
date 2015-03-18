

## To run:         source("/Users/alexgw/R/src/agwFiguresAll.R")      ;

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

SOURCE_ALL_WHERE  <- "/Users/alexgw/R/src/"
FIGURES_ALL_WHERE <- "/Users/alexgw/R/"
sr <- function(file) {
     agwSrc(paste(SOURCE_ALL_WHERE, file, sep=''))
}


srmain <- function(file) {
     agwSrcAndRun(paste(SOURCE_ALL_WHERE, file, sep=''))
}


df <- function(folder) {
     dir.create(paste(FIGURES_ALL_WHERE, folder, sep=''), showWarnings=FALSE, recursive=TRUE, mode="0775")
}





#df("Figure2b") ; srmain("agwFigure2b.R") ;
#df("Figure3a") ; srmain("agwFigure3a.R") ;
#df("FigureA_Indicators")    ; srmain("agwFigureA_DetailPlots.R") ;
#df("FigureB_ModSummary")    ; srmain("agwFigureB_ModuleSummary.R") ;
#df("FigureC_Drug_Clusters") ; sr("agwFigureC.R") ;
df("Figure3b_OverallPosRecovery") ; srmain("agwFigure3b_OverallPos.R") ;
#df("FigureD_PredSummary") ; df("FigureD2_PredSummary") ; sr("agwFigureD_GeneralEval.R") ;



##
