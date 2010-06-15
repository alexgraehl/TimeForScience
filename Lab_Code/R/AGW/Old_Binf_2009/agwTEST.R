
##       source("/Users/alexgw/R/src/agwTEST.R");        ## <-- To reload this file

## Here is an example where the "drugs that target the module" set (the "x" variable / first argument in t.test) causes significantly more lethality to the particular ORF deletion than "drugs that target other modules" (the "y" variable / second argument in t.test)
##ttt <- t.test(c(15, 15, 15, 15, 15), c(12,11,12,13,14,12,15,16,12,11,11,8,12,12,11,10),alternative="greater")
## Just for testing

source("/Users/alexgw/R/src/agwStandardFunctions.R");
source("/Users/alexgw/R/src/agwDrugCalcFunctions.R");

if (length(dev.list()) > 0) { dev.off() }
#plot.new()

#agwGlobalLoad("GLOBAL.HOM.SENS.COMPENDIUM")
#agwGlobalLoad("GLOBAL.POSITIVES.COMPLETE.DATAFRAME")


zom <- agwNewHash()
agwHashOfHashesPut(zom, "z", "xxj", 123)

agwHashOfHashesPut(zom, "z", "y", 123)

print(agwHashOfHashesGet(zom, "z","y"))

print(agwHashOfHashesGet(zom, "z","x"))

print(agwHashOfHashesGet(zom, "z","xxj"))

#plot(hr)

#(GLOBAL.HOM.SENS.COMPENDIUM))

#a <- GLOBAL.HOM.SENS.COMPENDIUM["Hillenmeyer/hydroxyurea_50000_um#2"]
#a <- GLOBAL.HOM.SENS.COMPENDIUM["Hillenmeyer/hydroxyurea_50000_um#2"]
#b <- GLOBAL.HOM.SENS.COMPENDIUM["Hillenmeyer/benomyl_10_ug#1"]
#b <- GLOBAL.HOM.SENS.COMPENDIUM["Hillenmeyer/hydroxyurea_50000_um#1"]

#print(  cor(a,b, method="pearson", use="complete.obs" ) )

#cor(
#    GLOBAL.HOM.SENS.COMPENDIUM["Hillenmeyer/cisplatin_500_um#1"]

    ## [181] "Hillenmeyer/hydroxyurea_100000_um#1"
    ## [182] "Hillenmeyer/hydroxyurea_100000_um#2"
    ## [183] "Hillenmeyer/hydroxyurea_200000_um"
    ## [184] "Hillenmeyer/hydroxyurea_25000_um"
    ## [185] "Hillenmeyer/hydroxyurea_50000_um#1"
    ## [186] "Hillenmeyer/hydroxyurea_50000_um#2"
    ## [187] "Hillenmeyer/hydroxyurea_50000_um#3"
    ## [188] "Hillenmeyer/hygromycin_0.4_um"
    


##options(error=recover)
## (getenv "USER")
##rm(gv.HOM.RANKS.BY.DRUG)
##rm(gvCompleteCollection)
##rm(gvMemberStatsList)








  ##ss <- matrix(c(33,22, 11,33, 22,11, 55,55, 44,44),ncol=2,byrow=TRUE) ## <-- for testing how to sort matrices
  ##rr <- matrix(c("YES","YES", 0,0, 0,0, 0,0 , "YES","YES"),ncol=2,byrow=TRUE)
  ##apply(ss,BY.COL,function(a){ iord <- order(a) ;  rr[iord]  } )



                                        #agwSensVsOverlapPlot(moduleName=names(gvDrugTargetMapping)[1]                     , theStatCollection=gvNeighborStatsList                     , theModuleMemberCollection=gvMemberStatsList                     , filePath=agwGlobalOutPath("")                     , fileName="S_test") ; print("ENDING EARLY!!!") ; stop("ENDED EARLY--DEBUGGING") ##"OUR_KNOWN_TARGET:HAND_CURATED:phosphatase_inhibitor"
###


#  source("/Users/alexgw/R/src/agwFigure2b.R");

print("Finished loading and running agwTEST.R")
