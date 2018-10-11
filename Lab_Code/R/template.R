#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#source("~willia51/workspace/0_CODE/common/common.R")
#source("~willia51/workspace/0_CODE/common/multi_assay_experiment.R")
#source("~willia51/workspace/0_CODE/common/paths.R")
#source("~willia51/workspace/0_CODE/common/salmon_quant.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (interactive()) { options(error=recover, max.print=1000) } else { options(error=NULL) } # useful even non-interactively
# Error = traceback DOES NOT STOP THE SCRIPT
options(stringsAsFactors=F, menu.graphics=F) # for top of file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Common functions
print0<-function(...){print(paste0(...))};
system0<-function(...){m=paste0(...);print0("[SYSTEM]: ",m);system(m)}
mandate<-function(x,...){if(!x){msg=paste0(...);warning(msg);stop(msg)}}
file.nonzero.exists <-function(f){return(file.exists(f)&&file.info(f)$size>0)}
first_existing <- function(...){for(x in list(...)){if(file.exists(x)){return(x)}};return(NA)}
devclear       <- function() { while (!is.null(dev.list())) { dev.off() } }
GLOBAL_ERRORS <- c("")
errlog <- function(...) { msg=paste0(...);print0(msg);warning(msg); GLOBAL_ERRORS <<- append(GLOBAL_ERRORS, msg); }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr); library(tibble); library(readr);
library(MultiAssayExperiment); library(SummarizedExperiment); # Note that if you include MAE, you should also inclued SummarizedExperiment so that e.g. 'names(...)' works
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Convenient functions for making tables in knitr documents
# Remember to enclose in an .Rmd in ```{r results="markdown"} ... ```
agw_dtfile    <- function(file, caption=NULL, ...) {
     require("readr")
     if (is.null(caption)) { caption = paste0("Data table: ", basename(file)) }
     if (grepl("[.](tsv|tab)", file, ignore.case=T)) { x = readr::read_tsv(file) }
     else if (grepl("[.]csv", file, ignore.case=T))  { x = readr::read_csv(file)  }
     else { warning("Couldn't guess the file type from the extension (trying tab-delim)."); x = readr::read_tsv(file) }
     agw_datatable(data=x, caption=caption, ...)
}
agw_datatable <- function(data, caption="Data Table", pageLength=10, options=NULL, ...) {
  require("DT"); # install.packages("DT")
  if (is.null(data) || is.na(data)) { stop("Your data to 'agw_datatable' was NULL or NA!") }
  if (is.null(options)) { options=list(pageLength=pageLength, lengthMenu=c(pageLength,50,200,999), autoWidth=TRUE) }
  numeric_colnames = colnames(data)[sapply(data, function(ccc) { return(is.numeric(ccc) && !is.integer(ccc))})] # find which columns are NUMERIC but not integers so we can round them
  SIGNIF_DIGITS = 4
  DT::datatable(data=data, caption=caption, filter=list(position='top',plain=TRUE),  rownames=FALSE, options=options, ...) %>% formatSignif(numeric_colnames, SIGNIF_DIGITS)   #style="bootstrap", . Note: works fine with NULL (c()) numeric_colnames input
}

args        <- base::commandArgs(trailingOnly=TRUE)
if (length(args) > 0 && !grepl("^[-]", args[1])) {
     stop(paste0("[:ERR:] Usage error: It looks like you are using this script by passing in something BESIDES an option, which must start with a '-'. You can't do that: try running this script with '--help' for usage information."))
}
library(optparse)
option_list <- list(optparse::make_option(c("-i", "--input"), dest="input", type="character", help="Input SAMPLE_LIST config file, which is actually an R list of lists!", metavar="FILE"),
                    optparse::make_option(c("-s", "--something"), dest="somethingvar", type="character", help="Another sample", metavar="SOMETHING")
                    )
opt        <- optparse::parse_args( optparse::OptionParser(option_list=option_list) ) #opt_parser)

if (!is.null(opt$input)) {
     print(opt$input)
} else {
     print("Apparently the '--input' was not specified on the command line.")
}

# Don't load all this stuff until AFTER you handled the argument parsing! It's slow.
#library(reshape2)
#library(pheatmap)
# library(ggplot2)
#library(viridis); # viridis::magma() / viridis::inferno() / viridis() color schemes. install.packages('virids')
library(readr); library(tibble); library(dplyr)
library(pryr); # pryr::object_size(obj1, obj2, ...)

print("[Done]")
