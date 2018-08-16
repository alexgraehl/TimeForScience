
# ======================================================================================
# First, attempt to load Alex's common utility code off the filesystem. Hopefully this already exists somewhere!
# Note that we NORMALLY attempt to load this based on where the "R_BINF_CORE" variable is, but we also a few hard-coded locations in case that environment variable isn't there.
for (tryFile in c(file.path(Sys.getenv("R_BINF_CORE"), "Utility", "agwUtil.R")
                , "/work/Common/Code/R_Binf_Core/Utility/agwUtil.R"
                , "~/TimeForScience/Lab_Code/R/AGW/agwUtil.R", "./agwUtil.R")) {
     if (file.exists(tryFile)) { source(tryFile); break; }
}
if (!exists("print.agw")) {
     if (0 == nchar(Sys.getenv("R_BINF_CORE"))) {
          warning("----------------------------------------\nERROR: The environment variable R_BINF_CORE was not defined.\nThis variable must be set to wherever the R Bioinformatics Core scripts root directory is.\nYou can set R_BINF_CORE by adding the following line to the end of your ~/.bashrc file:\n      export R_BINF_CORE=/SOME/PLACE/WHERE/R_Binformatics_Core_Libs_are_located <-- obviously not the real path\nThen you can open a new shell and run this script.\n\nIf you are running this somewhere without the shell set up, and want to hard-code the location in the file, you can also set it on the command line with: Sys.setenv(R_BINF_CORE=\"/work/Common/Code/R_Binf_Core\")\n----------------------------------------\n")
     }
     warning("Error -- can't find Alex's R utilities on the filesystem.")
}

# ======================================================================================
if (interactive()) { options(error=recover, max.print=1000) } else { options(error=NULL) } # useful even non-interactively
# Error = traceback DOES NOT STOP THE SCRIPT
options(stringsAsFactors=F, menu.graphics=F) # for top of file
# ======================================================================================
# Common functions
print0<-function(...){print(paste0(...))};
system0<-function(...){m=paste0(...);print0("[SYSTEM]: ",m);system(m)}
mandate<-function(x,...){if(!x){msg=paste0(...);warning(msg);stop(msg)}}
file.nonzero.exists <-function(f){return(file.exists(f)&&file.info(f)$size>0)}
first_existing <- function(...){for(x in list(...)){if(file.exists(x)){return(x)}};return(NA)}
devclear       <- function() { while (!is.null(dev.list())) { dev.off() } }
GLOBAL_ERRORS <- c("")
errlog <- function(...) { msg=paste0(...);print0(msg);warning(msg); GLOBAL_ERRORS <<- append(GLOBAL_ERRORS, msg); }
library("dplyr"); library("readr"); library("tibble")
# ======================================================================================


