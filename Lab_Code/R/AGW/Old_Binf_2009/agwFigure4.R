##       source("/Users/alexgw/R/src/agwFigure4.R");        ## <-- To reload this file

k_FIGURE_4_FOLDER_NAME <- "Figure4"

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");


## =========================================================================

## =========================================================================

agwFigure4Plot <- function() {
     agwGlobalLoad("gv.HOM.RANKS.BY.DRUG")
     agwGlobalLoad("gv.HOM.RANKS.BY.ORF")
     agwGlobalLoad("gvNeighborStatsHash")
     pMat <- matrix() # gene vs pathway
     mMat <- matrix() # gene vs pathway membership
     
     drugNamesVec <- agwHashKeys(gvNeighborStatsHash)

     for (d in drugNamesVec) {
          thisDrugHash <- agwHashGet(gvNeighborStatsHash, d)
          orfNamesVec  <- agwHashKeys(thisDrugHash)
          
     }
}

agwFigure4_main <- function() {
     print("=========================================")
     print("About to plot Figure 4A for all drugs...")
     print(agwGlue("Going to save in the Figure 4A directory named <", k_FIGURE_4_FOLDER_NAME, ">..."))
     
     agwFigure3aPlot()
}
