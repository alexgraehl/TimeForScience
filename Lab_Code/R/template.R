
# ======================================================================================
# First, attempt to load Alex's common utility code off the filesystem. Hopefully this already exists somewhere!
# Note that we NORMALLY attempt to load this based on where the "R_BINF_CORE" variable is, but we also a few hard-coded locations in case that environment variable isn't there.
for (tryFile in c(file.path(Sys.getenv("R_BINF_CORE"), "Utility", "agwUtil.R")
                  , "/work/Common/Code/R_Binf_Core/Utility/agwUtil.R"
                  , "~/TimeForScience/Lab_Code/R/AGW/agwUtil.R"
                  , "./agwUtil.R")) {
     if (file.exists(tryFile)) { source(tryFile); break; }
}
if (!exists("print.agw")) {
     if (0 == nchar(Sys.getenv("R_BINF_CORE"))) {
          warning("----------------------------------------\nERROR: The environment variable R_BINF_CORE was not defined.\nThis variable must be set to wherever the R Bioinformatics Core scripts root directory is.\nYou can set R_BINF_CORE by adding the following line to the end of your ~/.bashrc file:\n      export R_BINF_CORE=/SOME/PLACE/WHERE/R_Binformatics_Core_Libs_are_located <-- obviously not the real path\nThen you can open a new shell and run this script.\n\nIf you are running this somewhere without the shell set up, and want to hard-code the location in the file, you can also set it on the command line with: Sys.setenv(R_BINF_CORE=\"/work/Common/Code/R_Binf_Core\")\n----------------------------------------\n")
     }
     warning("Error -- can't find Alex's R utilities on the filesystem.")
}

# ======================================================================================
options(stringsAsFactors=F, menu.graphics=F, error=recover) # for top of file
# ======================================================================================
# Common functions
print0=function(...){print(paste0(...))}; system0=function(...){print0("[SYSTEM CALL]: ",...);system(paste0(...))}
file.nonzero.exists=function(f){return(file.exists(f)&&file.info(f)$size>0)}
GLOBAL_ERRORS <- c("")
errlog <- function(...) { msg=paste0(...);print0(msg);warning(msg); GLOBAL_ERRORS <<- append(GLOBAL_ERRORS, msg); }
qsubize <- function(cmd, name_prefix="X") { # Returns a qsub wrapper for the command 'cmd'
     QSUB_EXE <- system2("which", args=c("qsub"),  stdout=T); stopifnot(file.exists(QSUB_EXE))
     username <- Sys.info()[["user"]]
     persistent_n <- ifelse(is.null(attr(qsubize, "sum")), yes=0, no=(attr(qsubize, "sum")+1)) # <-- 'n' will be a "static" variable that persists between function calls
     attr(qsubize, "sum") <<- persistent_n # note the "<<-" to save this to the global scope!!!
     cmd <- paste0("echo \"", cmd, "\" | ", QSUB_EXE, " -V -N \"B2B_", name_prefix, "_", persistent_n, "_", username, "\"", "\n")
     return(cmd)
}
# ======================================================================================

