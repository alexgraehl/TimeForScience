
## This file should not *source* any others!
## It is COMPLETELY SELF CONTAINED. So you could copy and paste all this code anywhere else.

## ==============================================
### Recommended header code template for other files to include this file ###
#for (tryFile in c(file.path(Sys.getenv("R_BINF_CORE"), "Utility", "agwUtil.R"), "/work/Common/Code/R_Binf_Core/Utility/agwUtil.R", "~/TimeForScience/Lab_Code/R/AGW/agwUtil.R", "./agwUtil.R")) {
#     if (file.exists(tryFile)) { source(tryFile); break; }
#}
#if (!exists("print.agw")) {
#	if (0 == nchar(Sys.getenv("R_BINF_CORE"))) {
#	     warning("----------------------------------------\nERROR: The environment variable R_BINF_CORE was not defined.\nThis variable must be set to wherever the R Bioinformatics Core scripts root directory is.\nYou can set R_BINF_CORE by adding the following line to the end of your ~/.bashrc file:\n      export R_BINF_CORE=/SOME/PLACE/WHERE/R_Binformatics_Core_Libs_are_located <-- obviously not the real path\nThen you can open a new shell and run this script.\n\nIf you are running this somewhere without the shell set up, and want to hard-code the location in the file, you can also set it on the command line with: Sys.setenv(R_BINF_CORE=\"/work/Common/Code/R_Binf_Core\")\n----------------------------------------\n")
#	     stopifnot(nchar(Sys.getenv("R_BINF_CORE")) > 0) ## <-- Make sure the R_BINF_CORE environment variable is defined (should be in your .bashrc file as export R_BINF_CORE=/some/path/to/Binformatics_Core_R/ ) and non-blank!
#	}	
#	stop("Error -- can't find Alex's R utilities on the filesystem.")
#}
#options(stringsAsFactors=FALSE)
#options(error=recover)

## ================ END OF HEADER TEMPLATE ==============================

cols <- Sys.getenv("COLUMNS") ## Get UNIX terminal column width...
if (nzchar(cols)) { options(width = as.integer(cols)) } # set the current display column width to whatever the number of columns in the terminal is!

Z_LINE <- '--------------------------------------------------------------------------------\n' ## 80 characters

## ====================================
# "Constants" for mtext and other plotting commands
SIDE_BOTTOM <- 1 ; SIDE_LEFT   <- 2 ; SIDE_TOP    <- 3 ; SIDE_RIGHT  <- 4
## ====================================

LAS_AXIS_PARALLEL      <- 0 ; LAS_AXIS_HORIZONTAL    <- 1
LAS_AXIS_PERPENDICULAR <- 2 ; LAS_AXIS_VERTICAL      <- 3
## ====================================

APPLY_BY_ROW <- 1 ; APPLY_BY_COL <- 2 ; APPLY_BY_BOTH <- c(1,2)
BY.ROW = 1; BY.COL = 2 ; BY.BOTH = c(1,2)
## ====================================

PCH.X = 4
PCH.DIAMOND = 23  ## The hollow diamong plotting character

PCH.SOLID.BOX = 15
PCH.BOX       = PCH.SOLID.BOX

PCH.HOLLOW.UP.TRIANGLE    = 2
PCH.HOLLOW.DOWN.TRIANGLE  = 6
PCH.BORDER.UP.TRIANGLE    = 24
PCH.BORDER.DOWN.TRIANGLE  = 25


PCH.SMALL.CIRCLE  = 20
PCH.MEDIUM.CIRCLE = 16
PCH.BIG.CIRCLE    = 19

PCH.SOLID.BOX = 15
PCH.BOX = PCH.SOLID.BOX

PCH.FILL.BOX           = 15
PCH.FILL.MEDIUM.CIRCLE = 16
PCH.FILL.UP.TRIANGLE   = 17
PCH.FILL.DIAMOND       = 18
PCH.FILL.BIG.CIRCLE    = 19
PCH.FILL.SMALL.CIRCLE  = 20

PCH.COLOR.CIRCLE  = 21
PCH.COLOR.SQUARE  = 22
PCH.COLOR.DIAMOND = 23
PCH.COLOR.UP.TRIANGLE   = 24
PCH.COLOR.DOWN.TRIANGLE = 25

LTY.DASHED.LINE     = 3  ## "LTY: Line TYpe"
LAS.HORIZONTAL.TEXT = 1  ## Text alignment for margin text
LAS.PERPENDICULAR   = 2
LAS.VERTICAL.TEXT   = 3  ## Text alignment for margin text

POS.BOTTOM = 1; POS.LEFT = 2; POS.ABOVE = 3 ; POS.RIGHT = 4 ;

A.X.AXIS = 1
A.Y.AXIS = 2

WORST.P.VALUE = 1

HIGH.RES.DPI       <- 150
SCREEN.RES.DPI     <- 72

DEFAULT.POINTSIZE  = 12

STANDARD.PLOT.WIDTH  = 750
STANDARD.PLOT.HEIGHT = 750

MAKE.AGW.SUBDIR = "make.tmp"

## ====================================
join.left.outer.agw <- function(left.mat, right.mat, both.key, left.key, right.key, silent=FALSE) {
     ## This is like an SQL "left outer join" . Look it up on google for more details!
     ## It always returns the data from left.mat, no matter whether it has a key in right.mat or not.
     ## Accepts the special name "rownames()", which means that instead of looking for a column index,
     ## we use the rownames

     ## This is pretty much just a re-invention of "merge" .
     
     SPECIAL_CASE_ROWNAME_STRING <- "rownames()"
     
     if (!missing(both.key)) { ## both.key was specified...
          assert.agw(missing(left.key) && missing(right.key), "join.left.outer.agw: You cannot specify both.key AND ALSO left.key or right.key. Both.key sets BOTH of the other keys!") ## both of the other key types must be missing...
          left.key <- right.key <- both.key
     }
     
     theVectorForMatching <- function(theMatrix, theKey) {
          ## If the key is the literal string "rownames()", then there is no column with the key,
          ## the key column is already the rownames!
          if (theKey == SPECIAL_CASE_ROWNAME_STRING) {    return(rownames(theMatrix))  }
          else                        {    return(theMatrix[, theKey])  }
     }
     
     if (!silent) { print.agw("Join in progress on keys <", left.key, "> and <", right.key, ">...") }
     matchingIDs <- match(  theVectorForMatching(left.mat,  left.key)
                          , theVectorForMatching(right.mat, right.key))
     
     if (right.key != SPECIAL_CASE_ROWNAME_STRING && left.key != SPECIAL_CASE_ROWNAME_STRING) {
          ## Neither key is "rownames", so let's not double-add this key
          colIndexToOmit <- which( colnames(rightMatch) == right.key)
          assert.agw(length(colIndexToOmit) == 1)
          ## don't double-add the same key to the resulting matrix...
          rightMatch <- right.mat[matchingIDs, -c(colIndexToOmit)]
     } else {
          rightMatch <- right.mat[matchingIDs, ]
     }
     
     if (!silent) { print.agw("Join operation completed.") }
     return(cbind(left.mat, rightMatch))
}



## ==============================================
## Makes a new directory, or a bunch of directories. (Input is a single name, or a vector.)
## Doesn't warn you if the directory exists, but does abort with an error
## if it fails to create the directory. (Unless ignore.errors=TRUE, in which case it continues anyway)
## ==============================================
mkdir.agw <- function(..., ignore.errors=FALSE) {
     for (d in c(...)) {
          if (is.na(d) || is.null(d) || nchar(d) == 0) {
               # Skipping creating a directory with invalid inputs...
          } else {
               dir.create(d, recursive=TRUE, showWarnings=FALSE)
               # showWarnings is FALSE above, because otherwise a warning is generated
               # every time the directory *already* exists, which is not a problem, and
               # does not warrant cluttering the screen with a superfluous warning.
               if (!file.exists(d)) {
                    warning(paste("mkdir.agw: Unable to create the new directory", d))
                    if (!ignore.errors) { stop("Failed to create a directory: halting.") }
                    # But if we fail to make the directory, that actually probably IS a problem.
               }
          }
     }
}


## ==============================================
## Prints a warning message and then exits, if the "thingThatShouldBeTrue" is not true.
## Much like C's assert
## ==============================================
assert.agw <- function(assertionAGW=NULL
                       , ...
                       , ignore.errors=FALSE
                       , use.browser=FALSE ## jump directly into the browser
                       , sep='') {

     if (length(assertionAGW) > 1) {
          warning("assert.agw: Assertions cannot be more than one element in length. Probably you need to wrap your assertion in the \"all(...)\" command. Perhaps you were checking to see if all elements of a list were equal.")
     }
     
     if (is.null(assertionAGW) || (length(assertionAGW) <= 0) || !assertionAGW) {
          if (!missing(...)) {
               print.red.agw(paste("Assertion failed: ", ..., sep=sep))
               warning(paste("Assertion failed: ", ..., sep=sep))
          } else {
               theStr <- ">>> Assertion failed, but no assertion error message was specified."
               print.red.agw(theStr)
               warning(theStr)
          }
          
          if (!ignore.errors) {
               if (use.browser) {
                    browser()
               } else {
                    stop("Assertion failed! Stopping now.")
               }
          } else {
               warning("continuing...")
          }
          
     }
}


## ==============================================
## t.test, but returns NA instead of stopping with an error if you give it values
## that don't make sense for a t.test (for example: variance of zero data)
## ==============================================
t.test.no.errors.agw <- function(...) {
     obj <- try(t.test(...), silent=TRUE)
     if (is(obj, "try-error")) { return(NA) } else { return(obj) }
}


## ==============================================
## From: https://stat.ethz.ch/pipermail/r-help/2008-February/154172.html
## ==============================================
try.or.return.default.agw <- function(tryThis, default=NA) {
     result <- default
     tryCatch(result <- tryThis
              , error=function(e) {}) # don't do anything with errors...
     return(result)
}

## ==============================================
## From: https://stat.ethz.ch/pipermail/r-help/2008-February/154172.html
## ==============================================
fail.with.default.agw <- function(default=NULL, theFunc, ...) {
     function(...) try_default(theFunc(...), default)
}


## ==============================================
## Prints a message, the arguments to it pasted together (with no spaces between items by default)
## If you say "silent=TRUE" then it does nothing---this can be used to control printing for debugging
## ==============================================
print.agw <- function(...
                      , sep=''
                      , silent=FALSE
                      , newline=TRUE
                      , log=FALSE) {
     if (!silent) {
          ####### Print directly to the screen
          cat(paste(date(), "\t", sep=''))
          cat(..., sep='')           ## Print the main part of the data we wnat to print
          if (newline) { cat("\n") } ## Print a newline (maybe)
          ####### Done printing directly to the screen
          
          if (!is.null(log) && log != FALSE) {
               #### Write to a log file
               ## If "log" is true, then we also log to a file, as well
               ## as printing to stdout. If log is TRUE, then we print to
               ## a file named "agw.log.txt" . If log is a name, then we make
               ## a file with that name.
               if (log == TRUE) { fileToWriteTo = "agw.log.txt"; } ## Default name for the log file
               else { fileToWriteTo <- log } ## If "log" is something other than true or false, we assume it's the file name to write to
               cat(paste(date(), "\t", sep=''), file=fileToWriteTo, append=TRUE) ## <-- write to the log file!
               cat(..., sep='', file=fileToWriteTo, append=TRUE) ## <-- write to the log file!
               if (newline) { cat("\n", file=fileToWriteTo, append=TRUE) } # <-- Write a newline to the file
               #### Done writing to a log file
          }
          
     }
}


work.log.agw <- function(...) {
     print.agw(..., log="work.log.txt") ## A hard-coded "work log" file that would be suitable for researchers to read.
     print.agw(..., log=TRUE) ## Also log to the screen and to the regular log
}

log.agw <- function(..., log=TRUE) {
     print.agw(..., log=log)
}

## ==============================================
## Like "length" and "dim" together at last!
## Gives you some basic data on whatever variable you passed in.
## Also useful: "str"
## ==============================================
huh <- function(x) {
     if (is.null(x)) {
          print.magenta.agw("  The variable is NULL!")
     } else if (is.list(x)) {
          listSpacer <- "     "
          dataFrameStatus.text <- ifelse(is.data.frame(x), "(also a data frame) ", "")
          print.cyan.agw("  List with class ", paste(class(x), collapse=", "))
          print.cyan.agw("  List with attributes ", paste(attributes(x), collapse=", "))
          print.cyan.agw("  List ", dataFrameStatus.text, "of length ", length(x), ".\n")
          print.cyan.agw("Names in the list are:")
          print.cyan.agw(paste(listSpacer, paste(names(x), collapse=paste("\n", listSpacer, sep='')), sep=''))
     } else if (is.vector(x)) {
          print.green.agw("  Vector of type ", typeof(x))
          print.green.agw("  Length: ", length(x))
     } else if (is.matrix(x)) {
          print.yellow.agw("  Matrix of type ", typeof(x))
          print.yellow.agw("  Dim: ", nrow(x), " rows by ", ncol(x), " cols.")
          if (ncol(x) <= 50) { print.yellow.agw("  ", ncol(x), " cols: ", paste(colnames(x), collapse=", "))
          } else {             print.yellow.agw("  ", ncol(x), " cols (showing first 50): ", paste(c(colnames(x)[1:50], "..."), collapse=", ")) }
          print.yellow.agw("--------------------------------")
          if (nrow(x) <= 50) { print.yellow.agw("  ", nrow(x), " rows: ", paste(rownames(x), collapse=", "))
          } else {             print.yellow.agw("  ", nrow(x), " rows (showing first 50): ", paste(c(rownames(x)[1:50], "..."), collapse=", ")) }
     } else if (is.object(x)) {
          print.green.agw("  Object of class ", class(x)[[1]])
          print.green.agw("  With attributes ", paste(attributes(x), collapse=", "))
          print.green.agw("  With ", length(slotNames(x)), " slotNames(): ", paste(slotNames(x), collapse=', '))
          if (showMethods(class(x), printTo=FALSE)[3] == " <not a generic function>") {
               print.green.agw("  And with no methods.")
          } else {
               print.green.agw("  And with showMethods(): ", showMethods(class(x)))
          }
     } else {
          print.magenta.agw("  Type: ", typeof(x))
     }
}


## ==============================================
## Extracts a vector out of a list of lists. Loops over the OUTER list, and returns a particular element
## from each sub-list. For example, if you have a list named "ZZZ" of elements like
## list("a"=something, "b"=something),
## then you can get all of the "a"s out as a vector with listExtract.agw(ZZZ, "a")
## If you say "silent=TRUE" then it does nothing---this can be used to control printing for debugging
## ==============================================
listExtract.agw <- function(listOfLists, key) {
     ## maybe this is the same as "unlist?"
     return(sapply(listOfLists, "[[", key))
}
#sapply(a2,"[[","ensembl_ID")



## ==============================================
## Emulates the "make" utility.
## Only computes something once---on subsequent runs, the variable is actually
## loaded from a file instead of being recomputed.
               ## Does NOT save to the global namespace.
               ## in order to use this function, you MUST
               ## assign the return value.
## Example:  thing <- make.agw("myvar", makeTheThing(x=10) )
## This will compute thing using the "makeTheThing" function once, and save the result in
## "myvar.RData.tmp". Later it will just assign the "thing" variable from that file,
## rather than recomputing it.
## ==============================================
make.agw <- function(varname, thingToRunIfVariableIsNotFound, varToCheckInsteadOfVarname) {
     mkdir.agw(MAKE.AGW.SUBDIR)
     filename <- file.path(MAKE.AGW.SUBDIR, paste("tmp.make.", varname, ".RData.tmp", sep=''))
     
     if (!exists("z.agw.make.files")) { z.agw.make.files <<- c() }
     if (!exists("z.agw.make.vars")) {  z.agw.make.vars  <<- c() }
     z.agw.make.files <<- union(z.agw.make.files, filename)
     z.agw.make.vars  <<- union(z.agw.make.vars, varname)
     
     if (!file.exists(filename)) {
          print.green.agw("make.agw: Making the file <", filename, "> and saving the R variable <", varname, "> to it for later retrieval without having to recompute it.")
          ## If the file "filename" does NOT already exist,
          ## then we will run the "thingToRunIfVariableIsNotFound" function,
          ## and save the result into a file for future reading (without having to recompute it)
          temp <- thingToRunIfVariableIsNotFound
          save(temp, file=filename)
          return(temp)
     } else {

          variableAlreadyInNamespace <- FALSE
          
          if (!missing(varToCheckInsteadOfVarname)) {
               variableAlreadyInNamespace <- (exists("varToCheckInsteadOfVarname") && !is.null(varToCheckInsteadOfVarname) && !is.na(varToCheckInsteadOfVarname))
          } else {
               variableAlreadyInNamespace <- agwHasContent(varname)
          }
          
          if (variableAlreadyInNamespace) {
               ## If the variable already exists and isn't null / NA, then load it again...
               if (!missing(varToCheckInsteadOfVarname)) {
                    print.yellow.agw("make.agw: Skipping reloading of the variable we were supposed to check. varname was <", varname, ">,\nbut the variable to actually check was overriden in the call to varToCheckInsteadOfVarname as <", paste(substitute(varToCheckInsteadOfVarname), collapse=','), ">.\nYou can remake this variable by deleting the on-disk file at <", filename, ">.\n")
                    return(varToCheckInsteadOfVarname) ## <-- gotta return it, because the user is going to try to assign another variable based on this!
               } else {
                    print.yellow.agw("make.agw: Skipping reloading of variable <", varname, ">, as it already exists in the namespace.\nIf you want to remake it, delete the on-disk file at <", filename, ">.\n")
                    return(get(varname)) ## <-- gotta return it, because the user is going to try to assign another variable based on this!
               }
          } else {
               print.yellow.agw("make.agw: Retrieving data directly from <", filename, ">, without doing any additional computation.\nMake sure you save the return value of this function, as the loaded variable does NOT become part of the namespace.\n")
               ## If the file "filename" *IS* found already,
               ## then load the variable we want out of it, and return it.
               local({
                    ## Does NOT save to the global namespace.
                    ## in order to use this function, you MUST
                    ## assign the return value.
                    loadedVariables.vec <- load(filename)
                    assert.agw(length(loadedVariables.vec) == 1, "Currently, we can only save or load ONE variable at a time using make.agw")
                    return(get(loadedVariables.vec[1]))
               })
          }
     }
}

make.list <- function() {
    print(z.agw.make.vars)
    #print(z.agw.make.files)
}

make.clean <- function(varnames.vec) {
     if (missing(varnames.vec)) { ## Clean EVERYTHING if nothing is specified
          ## clean ALL
          varnames.vec <- z.agw.make.vars
     }
     
     if (length(varnames.vec) == 0) {
          print("Nothing to clean.")
          return()
     }
     
     sapply(varnames.vec,
            function(v) {
                 file <- file.path(MAKE.AGW.SUBDIR, paste("tmp.make.", v, ".RData.tmp", sep=''))
                 print.agw("Removing variable <", v, "> and file <", file, ">...")
                 if (exists(v)) {  rm(v)  }
                 if (file.exists(file)) { file.remove(file) }
            })
     z.agw.make.vars  <<- setdiff(z.agw.make.vars, varnames.vec)
     filenames.vec <- file.path(MAKE.AGW.SUBDIR, paste("tmp.make.", varnames.vec, ".RData.tmp", sep=''))
     z.agw.make.files <<- setdiff(z.agw.make.files, filenames.vec)
}

## ==============================================
## Same as "file.path", but requires that the file EXISTS ALREADY. (Otherwise, it quits with an error message!)
## First arguments ("..."): the first arguments are the elements to construct a path from (identical to the built-in R function file.path)
## missing.message (optional argument): An additional warning message to print out in the event that the file does NOT exist.
## ==============================================
file.path.that.exists.agw <- function(..., missing.message="File not found.") {
     paths.vec <- file.path(...)
     if (!all(file.exists(paths.vec))) {
          print.red.agw(missing.message)
          print.red.agw("Could not find a required file! 'file.path.that.exists.agw' was called. This function requires that every file passed in must already exist!")
          print.red.agw("Number of files that we could not find: ", sum(!file.exists(paths.vec)))
          print.red.agw("Filenames that could not be found: ", paste(paths.vec[!file.exists(paths.vec)], collapse=", "))
          assert.agw(all(file.exists(paths.vec)), missing.message)
     }
     return(paths.vec) # might be more than one path
}


## ==============================================
## Equivalent to asking "file.info$size > 0" and "file.exists" for each file.
## Returns a TRUE/FALSE vector of the same length as the number of input arguments.
## Example:
##   file.has.contents.agw("/bin", "/fake/thing.tmp")
## Result:
###  c(TRUE, FALSE)
file.nonzero.agw <- function(...) {
     if (!all(nchar(...) > 0)) { print("ERROR: You passed a ZERO LENGTH argument into 'file.nonzero.agw'!"); stop() }
     if (!all(!is.null(...))) { print("ERROR: You passed a NULL argument into 'file.nonzero.agw'!"); stop() }
     return(sapply(list(...)
                   , function(x) { return(file.exists(x) && file.info(x)>0) }
                   ))
}


## ==============================================
is.32.bit <- function() { return (4 == .Machine$sizeof.pointer) } ## Tells you if your copy of R is 32-bit
is.64.bit <- function() { return (8 == .Machine$sizeof.pointer) } ## Tells you if your copy of R is 64-bit
## ==============================================
## Similar to "source"
## ==============================================
## Provides a boolean list of indices to be removed from a vector /data frame / whatever.
## Throws an error if you have specified nonexistent names, unless you say ignore.mismatch=TRUE.
## Example:
# m  <- matrix(rnorm(100), ncol=10) ; colnames(m) <- letters[1:10]
# nn <- m[, omitByName.agw(colnames(m), c("a","i"))]  ## <-- removes columns "a" and "i"
## ==============================================
toOmitNames.agw <- function(allNames, omitThese, ignore.mismatch=FALSE) {
     to.keep.bool.vec <- !(allNames %in% omitThese)
     if (!ignore.mismatch) {
          itemsInRemoveThatAreAlsoInNames.bool.vec <- (omitThese %in% allNames)
          if (!all(itemsInRemoveThatAreAlsoInNames.bool.vec)) {
               print.magenta.agw("=========================================")
               print.magenta.agw("===================")
               print.red.agw("Error: The function \"", (match.call())[[1]], "\" was told to remove a name that did *not* exist in the input.")

               print.red.agw("       In other words, one of the items in omitThese was NOT also present in allNames.")
               print.red.agw("       Here is allNames:") ;     print(allNames)
               print.red.agw("       And here is omitThese:"); print(omitThese)
               print.red.agw("       This boolean vector shows which items in omitThese were found in allNames. The \"FALSE\" ones are bad news:") ; print(itemsInRemoveThatAreAlsoInNames.bool.vec)
               print.red.agw("       Thus, the items in omitNames that were NOT in allNames were:") ; print(omitThese[!itemsInRemoveThatAreAlsoInNames.bool.vec])
               print.magenta.agw("===================")
               print.magenta.agw("=========================================")
               assert.agw(all(itemsInRemoveThatAreAlsoInNames.bool.vec), "One of the names that we wanted to remove was *not* a name that was actually present. You may want to check your spelling AND capitalization, or override this with ignore.mismatch=TRUE.")
          }
     }
     return(to.keep.bool.vec)
}

## ==============================================
## Omits columns in a matrix or data frame by NAME.
## Throws an error if you have specified nonexistent names.
## (Remove columns by name instead of by index.)
## ==============================================
omitColumnsByName.agw <- function(mat, namesToRemoveVec, ignore.mismatch=FALSE) {
     return(mat[ , toOmitNames.agw(colnames(mat), namesToRemoveVec, ignore.mismatch=ignore.mismatch)])
}     

## =================================================================
## Shows the detailed information of all the sub-items of a list.
## It's like a more thorough version of "names(...)"
## =================================================================
agwNames <- function(thingToExamineInDetail) {
     v <- vector(mode="character", length=length(thingToExamineInDetail))
     if(typeof(thingToExamineInDetail) == "environment") {
          ## it's a hash/environment
          for (i in 1:length(thingToExamineInDetail)) {
               key = ls(thingToExamineInDetail)[i]
               valueType = typeof(thingToExamineInDetail[[key]])
               v[i] <- paste("Hash key #", i, "/", length(thingToExamineInDetail)," (", valueType, "): ", key, sep='')
          }
     } else {
          ## not a hash...
          for (i in 1:length(thingToExamineInDetail)) {
               thing = thingToExamineInDetail[i]
               thingName = names(thing)
               valueType = typeof(thing)
               v[i] <- paste("Item #", i, "/", length(thingToExamineInDetail)," (", valueType, "): ", thingName, sep='')
          }
     }
     return(v)
}


## =================================================================
## Calculate numbers that correspond to a particular character string. Sums together all the ascii (maybe ascii, anyway) numerical values for the letters in a string, and returns a final number. Then this number gets used by "fullColorGroups" to color the bars on the graph.
## Note that the "gsub('[\\d.]'..." part replaces any numbers or decimal points with a blank spot. This way various concentrations / replicates will normally have the same resulting string->sum value. However, this also causes some drugs with different names to map to the same color when their numeric parts are deleted. Oh well!
## Saturation goes from 0 to 1 (0 = grayscale)
## Value goes from 0 to 1 (0 = black, 1 = bright colors)
## =================================================================
agwColorsFromStrings <- function(inputVectorOfStrings, saturation=1, value=1, allowNumeric=FALSE) {
     substString <- '[\\d.]'
     if (allowNumeric) { substString <- '' }
     
     stringAsciiValueSums <- sapply(inputVectorOfStrings, function(sss) {
          sum(sapply(strsplit(paste("z",gsub(substString,'',sss,perl=TRUE)),'')[[1]], ## "paste("z") is only there to keep this from ever being empty!
                     function(x){ as.integer(charToRaw(x))}))
     })
     
     randuNCol = ncol(randu)  ## Randu is a built-in set of 400*3 uniform random numbers
     randuNRow = nrow(randu)
     fullColorGroups <- sapply(stringAsciiValueSums,function(a){ hsv(randu[1+((39*a) %% randuNRow),1+((17*a) %% randuNCol)], saturation, value) })
     return(fullColorGroups);
}
## =================================================================

## =================================================================
## Returns a multi-plot layout for however many items you want
## =================================================================
agwRectangularLayoutForNItems <- function(n) {
     stopifnot(length(n) == 1)
     ## Note! rows, cols. So it's the OPPOSITE of (x,y) <-- NOT x,y!!
     if (n <= 3) { return(c(n, 1))
     } else {
          NCOL <- ceiling(sqrt(n))
          NROW <- (1 * ((n %% NCOL) > 0)) + (n %/% NCOL)  # <-- %/% is "integer division" (x %/% y is like floor(x / y))
          return(c(NROW, NCOL))
     }
}

## =================================================================
## If data values cross the 0-axis, then return a ylim/xlim that will be symmetric
## around this axis.
## =================================================================
agwSymmetricAxis <- function(dataPoints) {
     top <- max(dataPoints, na.rm=T)
     bot <- min(dataPoints, na.rm=T)
     if (is.na(top)) { top=1; warning("No non-NA data points were passed to agwSymmetricAxis") }
     if (is.na(bot)) { bot=-1; warning("No non-NA data points were passed to agwSymmetricAxis") }
     if (top >= 0 && bot >= 0) { return(c(bot, top)) }
     if (top <= 0 && bot <= 0) { return(c(bot, top)) }
     else {
          absBiggest <- max(abs(top), abs(bot), na.rm=T)
          return(c(-absBiggest, absBiggest))
     }
}






## ====================================================
## ====================================================
## Prints with logging-to-file behavior by default
log.red.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="red", log=log)
}

log.yellow.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="yellow", log=log)
}

log.green.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="green", log=log)
}

log.cyan.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="cyan", log=log)
}

log.magenta.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="magenta", log=log)
}

log.blue.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="blue", log=log)
}
## ====================================================
## ====================================================

## ====================================================
## ====================================================
## Prints to the console, default is NOT to log to a file
print.red.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="red", log=log)
}

print.yellow.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="yellow", log=log)
}

print.green.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="green", log=log)
}

print.cyan.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="cyan", log=log)
}

print.magenta.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="magenta", log=log)
}

print.blue.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="blue", log=log)
}

## ====================================================
## ====================================================

print.color.agw <- function(..., newline=T, log=F, fg=NULL, bg=NULL) {
     # xterm256 is now deprecated. Use "xtermStyle" instead.
     if (require(xtermStyle)) {
          ## "require" checks the currently loaded packages and doesn't reload code that is already loaded.
          print.agw( xtermStyle::style(paste(..., collapse=NULL, sep=''), fg=fg, bg=bg) , newline=newline, log=log)
     } else {
          print.agw(..., newline=newline, log=log)
     }
}

## =================================================================

agwNames <- function(thingToExamineInDetail) {
     ## Shows the detailed information of all the sub-items of a list.
     ## It's like a more thorough version of "names(...)"
     v <- vector(mode="character", length=length(thingToExamineInDetail))
     if(typeof(thingToExamineInDetail) == "environment") {
          ## it's a hash/environment
          for (i in 1:length(thingToExamineInDetail)) {
               key = ls(thingToExamineInDetail)[i]
               valueType = typeof(thingToExamineInDetail[[key]])
               v[i] <- paste("Hash key #", i, "/", length(thingToExamineInDetail)," (", valueType, "): ", key, sep='')
          }
     } else {
          ## not a hash...
          for (i in 1:length(thingToExamineInDetail)) {
               thing = thingToExamineInDetail[i]
               thingName = names(thing)
               valueType = typeof(thing)
               v[i] <- paste("Item #", i, "/", length(thingToExamineInDetail)," (", valueType, "): ", thingName, sep='')
          }
     }
     return(v)
}


## =================================================================

agwGlue <- function(...) {
     return(paste(..., sep='')); ## Just "paste" with no separator
}
## =================================================================


agwPrint <- function(...) {
     cat(agwGlue(..., "\n"))
}

## =================================================================
## Clamps the values in a vector to certain minimums / maximums
agwClampVec <- function(x, min=NULL, max=NULL) {
     ## There is a probably a better way to do this...
     stopifnot(max <= min) # if (max < min) { stop("ERROR in arguments to agwClamp: max is greater than min!"); }
     if (!is.null(min)) { x[x < min] = min }
     if (!is.null(max)) { x[x > max] = max }
     return(x);
}
## =================================================================






## =================================================================
## Input: a tab-delimited file with a single row header AND a single column header.
## Set <row.names> to NULL to force auto-numbering of rows (useful if you don't care about the row names, and want to avoid the "duplicate 'row.names' are not allowed" problem)
read.file.into.data.frame.default.agw <- function(...) { return (agwReadFileIntoDataFrame(...)) } ## <-- old function name for the same thing
agwReadFileIntoDataFrame <- function(filename, allowRagged=FALSE, row.names=1, header=TRUE, sep="\t", comment.char="", check.names=FALSE, quote='') {
     return(read.table(filename, stringsAsFactors=FALSE
                       , header=header
                       , comment.char=comment.char
                       , sep=sep
                       , row.names=row.names
                       , quote=quote
                       , na.strings=c("NA","NaN","ND")
                       , check.names=check.names  ## Read row/col names verbatim, do not change them (by default) to be R variables!
                       , strip.white=TRUE
                       , allowEscapes=TRUE
                       , fill=allowRagged));
}
## =================================================================
## Reads from a file into a list of lists
## "row.name.column" lets you pick which column becomes the row name.
## Note that I think R requires that row names are unique.
## If "row.name.column" is set to NA, then no row names are set.
read.file.into.list.of.lists.agw <- function(...) { return(agwReadFileIntoListOfLists(...)) }  ## <-- old function name for the same thing
agwReadFileIntoListOfLists <- function(filename, sep="\t", row.name.column=NA) {
     initList <- readLines(filename)
     initList <- strsplit(initList, sep, fixed=TRUE)
     if (!is.na(row.name.column)) {
          ## Set the names to whatever is in the specified column (assuming of course that there is something there...
          ## Might not be unique? This could be a huge problem,
          ## so probably we should check for uniqueness of names
          names(initList) <- sapply(initList,
                                    function(a) {
                                         a[[row.name.column]]
                                    })
     }
     return(initList)
}

## =================================================================
read.file.into.matrix.agw <-  function(...) { return(agwReadFileIntoMatrix(...)) } ## <-- old function name for the same thing
agwReadFileIntoMatrix <- function(filename) {
     ## Input: a tab-delimited file with a single row header AND a single column header.
     ## Returns a matrix.
     return(as.matrix(read.table(filename, stringsAsFactors=FALSE, header=TRUE, comment.char="", sep="\t", row.names=1, quote="", na.strings=c("NA","NaN","ND"), check.names=FALSE)));
}
## =================================================================



## =================================================================

## Calculate numbers that correspond to a particular character string. Sums together all the ascii (maybe ascii, anyway) numerical values for the letters in a string, and returns a final number. Then this number gets used by "fullColorGroups" to color the bars on the graph.
## Note that the "gsub('[\\d.]'..." part replaces any numbers or decimal points with a blank spot. This way various concentrations / replicates will normally have the same resulting string->sum value. However, this also causes some drugs with different names to map to the same color when their numeric parts are deleted. Oh well!
## Saturation goes from 0 to 1 (0 = grayscale)
## Value goes from 0 to 1 (0 = black, 1 = bright colors)
agwColorsFromStrings <- function(inputVectorOfStrings, saturation=1, value=1, allowNumeric=FALSE) {
     substString <- '[\\d.]'
     if (allowNumeric) { substString <- '' }
     
     stringAsciiValueSums <- sapply(inputVectorOfStrings, function(sss) {
          sum(sapply(strsplit(paste("z",gsub(substString,'',sss,perl=TRUE)),'')[[1]], ## "paste("z") is only there to keep this from ever being empty!
                     function(x){ as.integer(charToRaw(x))}))
     })
     
     randuNCol = ncol(randu)  ## Randu is a built-in set of 400*3 uniform random numbers
     randuNRow = nrow(randu)
     fullColorGroups <- sapply(stringAsciiValueSums,function(a){ hsv(randu[1+((39*a) %% randuNRow),1+((17*a) %% randuNCol)], saturation, value) })
     return(fullColorGroups);
}
## =================================================================


## =================================================================

agwLine <- function(where, type) {
     ## Draws a two-colored line
     ## Example: agwLine(where=0.5, type='h') ## <-- horizontal line at 0.5
     ZERO.LINE.COLOR     = "blue"
     ZERO.LINE.COLOR.2   = "orange"
     
     if (type == 'v') {
          abline(v=c(where),col=ZERO.LINE.COLOR, lwd=3)
          abline(v=c(where), lty=LTY.DASHED.LINE,col=ZERO.LINE.COLOR.2, lwd=3)
     } else if (type == 'h') {
          abline(h=c(where),col=ZERO.LINE.COLOR, lwd=3)
          abline(h=c(where), lty=LTY.DASHED.LINE,col=ZERO.LINE.COLOR.2, lwd=3)
     } else {
          stop("Wrong type for agwLine");	
     }
}
## =================================================================


## =================================================================
agwNewHash <- function(size=100L) {
     ## Initialize a new hash
     ## You use it as follows:   myHash <- agwNewHash()
     return(new.env(hash=TRUE, parent=emptyenv(), size=size))
}
## =================================================================
## =================================================================
agwHashKeys <- function(hash) {
     ## Note: this is very slow for a large hash (100,000+ elements makes this take several seconds)
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     ls(env=hash) #### Returns a vector of all the hash keys
     ## To get the size of a hash, just use length(hash)
}
## =================================================================
## =================================================================
agwHashGet <- function(hash, key, notFoundValue=NULL) {
     ## Returns notFoundValue if the key was not found,
     ## otherwise returns the value in the hash
     ## This is constant-time, even with 500,000+ elements, unlike accessing
     ## vectors and lists by name.
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     if (is.null(key)) { return(notFoundValue) }
     stopifnot(length(key) == 1)
     if (nchar(key) == 0 || is.na(key)) { ## <-- blank key string -> return "not found"
          return(notFoundValue)
     } else {
          stopifnot(typeof(key) == "character")
          return(mget(key,hash,ifnotfound=list(notFoundValue))[[1]])
     }
}
## =================================================================

agwHashPut <- function(hash, key, value) {
     ### This is for setting a SINGLE key-value pair in a hash.
     ### The "getter" function is "agwHashGet" , which returns the value when given a key.
     ### Example:
     ### myHash <- agwNewHash()
     ### agwHashPut(myHash, "year", 1776)
     ### print(agwHashGet(myHash, "year"))
     ### Note: nothing happens if you try to add a "" / null key
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     if (length(key) == 0 || nchar(key) == 0) { return() }
     else { assign(key,value,hash) }
}
## =================================================================
agwHashOfHashesPut <- function(hash, key1, key2, value) {
     ## Puts it into a hash within a hash
     ## Example: agwHashOfHashesPut(zom, "z", "xxj", 123)
     stopifnot(typeof(hash) == "environment")
     subHash <- agwHashGet(hash, key1, notFoundValue=NULL)
     if (is.null(subHash)) {
          subHash <- agwNewHash()
          agwHashPut(hash=hash, key=key1, value=subHash)
     } else {
          stopifnot(typeof(subHash) == "environment")
     }
     agwHashPut(hash=subHash, key=key2, value=value)
}


agwHashOfHashesGet <- function(hash, key1, key2, notFoundValue=NULL) {
     ## Gets the value from a hash within a hash
     ## Example: print(agwHashOfHashesGet(zom, "z","y"))
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     subHash <- agwHashGet(hash, key1, notFoundValue=NULL)
     if (is.null(subHash)) {
          return(notFoundValue)
     } else {
          stopifnot(typeof(subHash) == "environment") ## check for a hash argument...
          return(agwHashGet(hash=subHash, key=key2, notFoundValue=notFoundValue))
     }
}


agwMatrixFromHashOfHashes <- function(hashOfHashes) {
     stopifnot(typeof(hashOfHashes) == "environment") ## check for a hash argument...
     # Take a hash of hashes and makes a matrix out of it.
     # The first set of keys becomes the... maybe rownames
     # the second set of keys is the other dimension's names.
     # The values become the contents.
     dim1KeyVec <- agwHashKeys(hashOfHashes)
     dim2KeyVec <- agwGetAllSecondTierHashKeys(hashOfHashes)
     
     mMat <- matrix(data=NA, nrow=length(dim1KeyVec), ncol=length(dim2KeyVec)) #, dimnames=list(dim1KeyVec, dim2KeyVec))
     browser()
     for (key1 in dim1KeyVec) {
          for (key2 in dim2KeyVec) {
               mMat[key1, key2] <- agwHashOfHashesGet(hashOfHashes, key1, key2, notFoundValue=NA)
          }
     }
     return(mMat)
}

# For a hash-of-hashes
agwGetAllSecondTierHashKeys <- function(hashOfHashes) {
     stopifnot(typeof(hashOfHashes) == "environment") ## check for a hash argument...
     # It's easy to get all the FIRST TIER keys in a hash of hashes.
     # But sometimes you want to get all the keys for the sub-hashes.
     # Returns a vector *set* (i.e., no duplicates) of all the keys
     dim1KeyVec <- agwHashKeys(hashOfHashes)
     dim2KeysHash <- agwNewHash()
     for (k in dim1KeyVec) {
          keysFor2 <- agwHashKeys(agwHashGet(hashOfHashes, k))
          agwHashPutMultiple(dim2KeysHash, keysFor2, keysFor2)
     }
     return(agwHashKeys(dim2KeysHash)) # <- returns a vector of any key that appeared as a key in the "second tier" of the hash of hashes
}

## =================================================================
## =================================================================
agwHashContains <- function(hash, key) {
     ### Returns FALSE if the key was not found, and TRUE otherwise.
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     if (nchar(key) == 0) { return(FALSE) }
     else { return(!is.null(mget(key,hash,ifnotfound=list(NULL))[[1]])) }
}
## =================================================================

agwHashPutMultiple <- function(hash, keys, values) {
     ### This is for adding MULTIPLE key-value pairs to a hash
     ### in one single call.
     ### Note: nothing happens if you try to add a "" / null key
     ### Also: note that you can add a VECTOR here as well if you want!
     ### myHash <- agwNewHash()
     ### agwHashPutPlural(myHash, c("year","food","color"), c(1776,"pizza","purple"))
     ### print(agwHashGet(myHash, "food"))
     ### print(agwHashGet(myHash, "color"))
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     if (length(keys) == 0) {
          return()
     } else {
          stopifnot(length(keys) == length(values))
          for (i in 1:length(keys)) {
               assign(keys[i],values[i],hash)
          }
     }
}


agwHashPutMembership <- function(hash, keys, setValue=TRUE) {
     ### Makes a hash that ONLY tells you if a key is a member or not.
     ### These hashes' keys do not have unique values--they are *all*
     ### set to whatever "setValue" is.
     stopifnot(typeof(hash) == "environment") ## check for a hash argument...
     if (length(keys) == 0) {
          return()
     } else {
          sapply(keys,
                 function(theKey) {
                      agwHashPut(hash=hash, key=theKey, value=setValue)
                 })
     }
}

## =================================================================

agwCollectionContains <- function(collection, element) {
     ### Note: O(n) checking through an (unsorted) array from the beginning for a match
     for (i in 1:length(collection)) {
          if (collection[i] == element) {
               #agwPrint(collection[i] , " was equal to ", element)
               return(TRUE)
          }
     }
     return(FALSE)
}


## =================================================================

mirrorYlim <- function(yValuesVec) { ## mirror / symmetrical / symmetry limits
     #### Input: the vector of Y values
     #### Output: an output vector c(min,max) where
     #### min is -max, min is < 0, and max is > 0.
     #### Usage: plot(x,y, ylim=mirrorYlim(y))
     low <- min(yValuesVec, na.rm=TRUE)
     high <- max(yValuesVec,na.rm=TRUE)
     if (!is.finite(low) || low >= 0) {  low <- -1 }
     if (!is.finite(high) || high <= 0) { high <- 1 }
     if (low > -high) { low  <- -high }
     if (high < -low) { high <- -low }
     return(c(low,high))
}

agwMirrorLim <- function(valuesVec) { ## mirror / symmetrical / symmetry limits
     return(mirrorYlim(valuesVec))
}
## =================================================================





## =================================================================
## =================================================================


agwMatrixWithoutMissingLines <- function(inputMatrix) {
     ## Given an input matrix, returns that same matrix minus any rows and columns
     ## that are *all* NA. As long as there is at least one non-NA element, we
     ## leave that row or column.
     return(inputMatrix[ which(apply(inputMatrix, BY.ROW, function(a) !all(is.na(a))))
                        ,which(apply(inputMatrix, BY.COL, function(b) !all(is.na(b))))])
}

## If you want to reduce a matrix in size by removing the rows
## and columns that don't have enough data points, then you should
## use this function.
## In the case of "minN=1", it is identical to "agwMatrixWithoutNALines"
agwMatrixLinesWithEnoughData <- function(inputMatrix, minN=1) {
     ## Given an input matrix, returns that same matrix minus any rows and columns
     ## without at least "minN" non-NA elements.
     return(inputMatrix[ which(apply(inputMatrix, BY.ROW, function(a) (sum(!is.na(a)) >= minN)))
                        ,which(apply(inputMatrix, BY.COL, function(b) (sum(!is.na(b)) >= minN)))])
}


agwMatrixFromLists <- function(inList) {
     ### Makes a vector with maybe-different-length lists.
     ### Pads the extra values with NA.
     ### All matrices must be of the same type
     ## Does NOT MATCH things up by row names!
     longestRowLength <- max(sapply(inList, length))
     tempVec <- vector()
     for (i in 1:length(inList)) {
          thisVecLen <- length(inList[[i]])
          v <- c(inList[[i]], rep(NA, longestRowLength-thisVecLen))
          tempVec <- append(tempVec, v)
     }
     return( matrix(tempVec, nrow=length(inList), ncol=longestRowLength, byrow=TRUE) )
}

## =================================================================
## Input: a LIST of vectors. So this would be valid: listOfVectors=list(c(a=1,b=2,c=3), c(a=4,z=5,f=6))
##        The input vectors have to have **unique names** for each element.
## Output: the above data would provide a merged matrix, like:
##       VEC1      VEC2
##  a     1         4
##  b     2        <NA>
##  c     3        <NA>
##  f    <NA>       6
##  z    <NA>       5
## Note that the order of rows is just the default sort order. Usually this means alphabetical sort.
## This is like "merge" but it accepts any number of arguments, and is a bit more limited in scope since it only does vectors.
agwMatrixMatch <- function(listOfVectors) {
     assert.agw(is.list(listOfVectors), "Input must be a list. The elements of this list must be vectors. i.e., list(c(a=1,b=2,c=3), c(a=4,z=5,f=6)) . Note that we use the names(...) of the vector to match things up, so these vectors should also all have non-null names.")
     assert.agw(all(sapply(listOfVectors, is.vector)), "All the items in the list of vectors have to be vectors!")
     assert.agw(all(sapply(listOfVectors, function(v) { !is.null(names(v)) })), "All the items in the list of vectors have to have names. But at least one of them did not have names associated with its values!")
     assert.agw(all(sapply(listOfVectors, function(v) { !any(duplicated(names(v))) })), "Dang, we don't know how to deal with DUPLICATED names in the vectors! All the names(...) must be unique.")
     allNames <- sort(unique(unlist(lapply(listOfVectors, names)))) ## every element in listOfVectors must be a vector! with names!
     theMat  <- matrix(NA, ncol=length(listOfVectors), nrow=length(allNames))
     rownames(theMat) <- allNames
     for (i in seq(from=1, to=ncol(theMat), by=1)) {
          theMat[, i] <- listOfVectors[[i]][allNames] ## note that any name requested that is NOT present will show up as 'NA' in the final 'theMat'
     }
     return(theMat)
}


agwMatrixMeanByCertainColumns <- function(in.mat, column.group.ids.vec) {
     ## groupmean / matrix by group mean / matrix average by column / groupaverage / colaverage
     ## in.mat: an input matrix
     ## column.group.ids.vec: a correspondence for each COLUMN to which group it's in.
     ## Example input:
     ## matrix:   GroupA  GroupA   GroupB   GroupC
     ##               1       3         7      NA
     ##               0       NA        8       8
     ## input column.group.ids.vec:
     ##    c("GroupA", "GroupA", "GroupB", "GroupC")
     ## output:
     ## matrix:   GroupA   GroupB   GroupC
     ##              2        7        NA
     ##              0        8         8        <-- replaced by MEANS, removing NA when required
     assert.agw(length(column.group.ids.vec == ncol(in.mat)), "Must be one group ID per column of the input matrix.")
     groupMeansPerGene.list <- list()
     uniqueGroupIDs <- unique(column.group.ids.vec)
     for (ggg in uniqueGroupIDs) {
          colIndicesThatMatchThisGroup <- which(column.group.ids.vec == ggg)
          assert.agw((length(colIndicesThatMatchThisGroup) > 0), "No groups matched? Weird, this should be impossible.")
          subsetThisGroupOnly.mat <- in.mat[ , which(column.group.ids.vec == ggg), drop=F]
          theseMeans.vec <- rowMeans(subsetThisGroupOnly.mat, na.rm=T)
          theseMeans.vec[!is.finite(theseMeans.vec)] <- NA
          groupMeansPerGene.list[[ggg]] <- theseMeans.vec
     }
     final.mat <- sapply(groupMeansPerGene.list, cbind)
     rownames(final.mat) <- rownames(in.mat)
     colnames(final.mat) <- uniqueGroupIDs
     return(final.mat) ## <-- now all the columns for the SAME group ID have been collapsed down to their mean value!
}



## =================================================================
agwPlotLinesViolinDensity <- function(inVec, col="#00000077",scaleToY=NULL, mirror=TRUE, border=NA) {
     ## Given an input vector of values, it calculates the
     ## density and plots a mirrored-across-y-axis density plot
     ## of the data points. "scaleToY" makes the density plot
     ## scaled such that the highest Y value is set to whatever
     ## "scaleToY" is.
     ##
     ## border = NA means "only draw the fill, not the lines"
     ## border = 1 means "draw a border"
     ## border = "red" means draw a border in the color red
     ## "lwd" controls the line thickness, but is not an option
     ## to this function.
     
     ## Must be appended to an existing plot!
     if (length(na.omit(inVec)) >= 2) {
          dens <- density(inVec)
          scaleFactor = 1
          if (!is.null(scaleToY)) {
               maxY = max(dens$y)
               scaleFactor = scaleToY / maxY
          }
          
          polygon(dens$x, dens$y * scaleFactor, col=col, border=border)
          if (mirror) {
               polygon(dens$x, -dens$y * scaleFactor, col=col, border=border)
          }
     } else {
          print("AGW Warning: Not enough data points to plot a density histogram.")
     }
}
## =================================================================

agwRemoveNAFromBoth <- function(xVec = NULL, yVec = NULL) {
     ## You pass in a list of x-coordinates and a list of y-coordinates.
     ## If *either* the x or the y is NA, then *both* the elements are skipped.
     ## Otherwise, both elements are included. In the end, a list consisting of the xVec and yVec
     ## where *neither* element was NA is returned.
     ## In other words:
     ## if x = c(1,2,NA, 3)
     ## and y = (c(NA, 4, 5, 6)
     ## Then this function will return a list(x=c(2,3),y=c(4,6))
     if (length(xVec) != length(yVec)) {
          agwPrint("First (X) Vector length: ", length(xVec))
          agwPrint("Second (Y) Vector length: ", length(yVec))
          stop("The vectors passed into agwRemoveNAForPlot(...) must be the same length!")
     }
     elementsToKeepVec <- (!is.na(xVec) & !is.na(yVec))
     return( list(  x = xVec[elementsToKeepVec]
                  , y = yVec[elementsToKeepVec] ) )
}

## =================================================================
## =================================================================

agwHasContent <- function(vName) {
     ## Pass in the NAME of a variable, not a variable!!
     ## Returns TRUE if there is some data that isn't NA or NULL
     ## in this structure or variable.
     ## Returns FALSE if the structure is length zero, if no
     ## such variable exists, or if the variable is set to NA or NULL.
     ## So here are some things returning false:
     ##  list()     NA     NULL     vector(length=0)
     stopifnot(typeof(vName) == "character")
     if (!exists(vName)) { return(FALSE); }
     v <- eval(parse(text=vName))
     if (length(v) == 0) { return(FALSE); }
     if (is.null(v) || is.na(v)) {     return(FALSE); }
     if (typeof(v) != "environment") { ## <-- hashes should not (and cannot) be checked for NA-ness
          if (length(v) == 1 && is.na(v)) { return(FALSE); }
     }
     return(TRUE);
}

## =================================================================


agwImageBreakLabelsFromColors <- function(whichColors, whichBreaks, digits=NULL) {
     ## Used for making a legend for the heatmap-style "image(...)"
     ## plot command.
     ## When you have a list of colors and breakpoints,
     ## and you want a legend, this is the function you call
     ## to get the labels for the legend.
     ## It gives you labels like "(1,3]" for a color, indicating
     ## that that color is used for values > 1 and <= 3.
     ## (i.e., legend("topleft", legend=agwImageBreakLabelsFromColors(...)...)
     theLabels     <- vector(length=length(whichColors))
     theLabels[1]  <- agwGlue("<= ", format(whichBreaks[2], digits=digits))
     for (i in 2:(length(whichColors)-1)) {
          theLabels[i] <- agwGlue("(", format(whichBreaks[i], digits=digits) , ", ", format(whichBreaks[i+1], digits=digits), "]")
     }
     theLabels[length(whichColors)] <- agwGlue(">", format(whichBreaks[length(whichColors)], digits=digits))
     return(theLabels)
}


## =================================================================

agwGetHeatColors <- function(n) {
     # Makes a heatmap-colored-like vector that I can change
     # later. I think it literally is basically heat.colors,
     # except that the reds on the low end are a bit more obviously
     # different.
     if (n == 1) { return(c("red")); }
     if (n == 2) { return(c("black","red")); }
     if (n == 3) { return(c("black","red","yellow")); }
     else { return(c("black", heat.colors(n)[c(1,4:n)], "white")) }
}



## =================================================================


## ==============================================
## system.agw: Basically equivalent to print(cmd) and then system(cmd)
## Inputs: a vector of items (or a single item) to be pasted together to make a single system command.
## So you can say: system.agw("ls ", FILES, " ; mkdir somedir ")
## If a command has four spaces in a row, that is treated as as "fake" linebreak for display purposes only.
## Note that items are concatenated WITHOUT any surrounding spaces, so you have to add YOUR OWN SPACES
## around arguments, so that they don't all run together.
## If you want to have spaces around your inputs by default, use sep=" ".
## By default, prints the command in FORMATTED fashion, with line breaks. If you don't like this,
## set formatted=FALSE.
## If "dryrun" is true, then it just PRINTS the command, but does not execute it.
## Returns the exit code of the command. Not sure what happens if there are compound commands!
## ==============================================
system.agw <- function(...
                       , sep=''
                       , formatted=TRUE # <-- if true, then we split commands with " | " and multiple spaces onto their own lines. Keeps the output cleaner, at the expense of potentially being confusing. Set to FALSE to just print the command verbatim.
                       , requireZeroExitCode=FALSE ## Aborts on a non-zero exit code
                       , dryrun=FALSE
                       , wait=TRUE
                       , log=FALSE ## Should we log to the output file?
                       , time=FALSE ## Should we print the time?
                       , bash=FALSE ## Should we use /bin/bash (FALSE means use the terrible default shell)
                       ) {
     
     theCommand <- paste(..., sep=sep) ## Default is to just mash all the input arguments together with no spaces.
     syscallPrefix <- "System Call: ";
     if (!wait)  { syscallPrefix <- "System (Background): " }
     if (dryrun) { syscallPrefix <- "System call dry run (not executed!): " }

     fgColor <- ifelse(dryrun, "#00FFFF" ## cyan = dry run
                       , "#00FF00") ## green = real run
     
     numPrefixSpaces <- 1 + nchar(syscallPrefix) # <-- how many spaces required to indent the following lines of the combined system call properly?
     prefixWhitespace <- paste(rep(' ', numPrefixSpaces), collapse='') ## <-- a whitespace region to indent things properly

     if (time) { print.color.agw(paste("Syscall started at [", Sys.time(), "] ", "\n", sep=''), fg=fgColor, log=log, newline=F) } ## <-- note the lack of a newline!
     
     if (formatted) { # "Formatted" means we should be auto-splitting the line onto new lines
          formattedCmd = gsub("    ", paste("\n", prefixWhitespace, sep=''), theCommand, perl=TRUE) ## <-- four spaces becomes a newline
          formattedCmd = gsub(" [|] ", paste("\n", prefixWhitespace,"| ", sep=''), formattedCmd, perl=TRUE) # <-- a vertical bar surrounded by spaces becomes a newline, then a vertical bar
          print.color.agw(syscallPrefix, formattedCmd, fg=fgColor, log=log)
     } else {
          print.color.agw(syscallPrefix, theCommand, sep='', fg=fgColor, log=log)
     }

     exitCode <- 0
     if (!dryrun) {
          if (is.null(theCommand) || (nchar(theCommand) == 0)) {
               warning("system.agw was called with a blank argument. We are not actually going to run this (nonexistent) command.)")
          }
          
          if (bash) {
               ## Run it through BASH instead of SH. This allows the '<( ... )' construction. (Example: >>># cat <(zcat somezippedfile.gz) #<<< <-- makes a "magical" temporary file to store the output of zcat somezippedfile.gz )
               exitCode <- system(paste("bash -c ", base::shQuote(theCommand, type='sh'), sep=''), wait=wait)  ## <-- the unix error/success code
          } else {
               ## Run it through the default shell, /bin/sh.
               exitCode <- system(theCommand, wait=wait)  ## <-- the unix error/success code
          }
     }

     if (time) { print.color.agw(paste("Syscall done at [", Sys.time(), "]", sep=''), fg=fgColor, log=log, newline=TRUE) }
     
     if (requireZeroExitCode) {
          assert.agw(exitCode == 0, paste("Uh oh! The exit code for the command was NOT zero (zero indicates a success in most utilities)! Instead, the exit code was: ", exitCode, "\n", sep=''))
     }
     
     return(exitCode)
}






### ===============================================================================
## Just a plain rectangle with a bunch of (wrapped) text in it.
## wraplen: the number of characters before we wrap text
## leftMartin: the left margin. Ranges from 0-1. 0 = no margin, 1 = figure is only margin, and no text can be printed
## topMargin: the top start location. Ranges from 0-1. 0 means "start at the top, no margin" and 1.0 means "don't show the text, it's all margin"
### ===============================================================================
textDescriptionPlotAGW <- function(text, wraplen=60, leftMargin=0.15, topMargin=0.01, cex=1) {
     prevMar <- par()$mar  ;  par(mar=c(0,0,0,0)) ## no margins for this plot...
     wordwrap <- function(str, len) {
          modStr <- (strsplit(str, "\n"))[[1]]
          return(paste(strwrap(modStr, width=len), collapse="\n"))
     }
     plot(c, axes=F, xlab=NA, ylab=NA, type='n', frame.plot=F)
     text(x=leftMargin, y=(1-topMargin), adj=c(0,1), cex=cex, labels=wordwrap(text, len=wraplen))
     par(mar=prevMar) ## restore the old margins...
}

plot.text.only.agw <- function(...) {
     textDescriptionPlotAGW(...)
}








## =================================================================
## Colors -- gradient, like heat.colors
## NOTE: YOU MAY WANT TO USE THE BUILT IN "colorpanel" function instead.
## =================================================================
colors.agw <- function(n=12, type="blueblackyellow2", reverse=FALSE) {
     # NOTE: YOU MAY WANT TO USE THE BUILT IN "colorpanel" function instead.
     # Returns a color gradient, much like "heat.colors(...)" Several options for "type" are available.
     stopifnot(is.numeric(n)); stopifnot(n >= 1)
     totalColors <- n
     colRange01 <- (0:(totalColors-1))/(totalColors-1)
     range1 <- (0:(ceiling(totalColors / 2)-1))/(ceiling(totalColors/2)-1) ## For two-color gradients (the left half of the colors)
     range2 <- (1:(floor(totalColors / 2)))/(floor(totalColors/2))         ## For two-color gradients (the right half of the colors)
     type <- tolower(type)
     if (type == "greenwhitered") { ## green -> white -> red. Don't use this is possible--use "blueblackyellow"! Colorblind people can't see this.
          col1 <- rev(hsv(h=rev(0.3+0.20*range1), s=range1, v=(1.0 - 0.2*range1)))
          col2 <- (hsv(h=rev(0.0+0.15*range2), s=range2, v=(1.0 - 0.15*range2)))
          col <- c(col1, col2)
     } else if (type == "greenblackred") { ## green -> black -> red . Don't use this is possible--use "blueblackyellow"! Colorblind people can't see this.
          col1 <- rev(hsv(h=rev(0.3+0.20*range1), s=rev(1.0-0.2*range1), v=range1))
          col2 <- (hsv(h=rev(0.0+0.10*range2), s=rev(1.0-0.2*range2), v=range2))
          col <- c(col1, col2)
     } else if (type == "blueblackyellow") { ## blue -> black -> yellow . Note: this is somewhat uglier than "blueblackyellow2"!
          col1 <- rev(hsv(h=rev(0.6+0.20*range1), s=rev(1.0-0.2*range1), v=range1))
          col2 <- (hsv(h=rev(0.15+0.10*range2), s=rev(1.0-0.2*range2), v=range2))
          col <- c(col1, col2)
     } else if (type == "blueblackyellow2") { ## blue -> black -> yellow, hand-picked. Looks better than blueblackyellow!
          left <- range1[2:length(range1)]
          col1 <- rev(hsv(h=0.65-0.1*left, s=1.0, v=left)) ## 0.65 to 0.55 = blue
          col2 <- (hsv(h=0.166, s=1.0, v=range2)) ## 0.166 = yellow
          col <- c(col1, "black", col2)
          #col <- c("#00C8FF", "#00AAF5", "#0082D7", "#2464A8", "#004064", "black", "#646400","#919114","#B6B61E","#D7D728","#FFFF00")
     } else if (grepl("^gr[ea]y", type)) { col <- gray(colRange01) }
     else if (grepl("^brown", type))   { col <- rev(hsv(h=rev(0.2*colRange01), s=colRange01, v=(1.0 - 0.7*colRange01))) }
     else if (grepl("^sepia", type))   { col <- rev(hsv(h=rev(0.3*colRange01), s=colRange01, v=rev(colRange01))) }
     else if (grepl("^heat", type))    {  ## Nicer heatmap colors
          if (totalColors == 2) {
               col <- c("#000066", "#FFFF99")
          } else if (totalColors >= 10) {
               col <- c("black", "#330033", "#440044", "#550055", "#770044", "darkred", heat.colors(n-6))
          } else {
               col <- c("#330033", "#770044", heat.colors(n-2))
          }
     } else if (grepl("^oldheat", type)) {
          col <- heat.colors(totalColors) ## R's very-bright heatmap colors
     } else {
          print(paste("An unrecognized color gradient type (<", type, ">) was passed into colorGradientAGW:", sep=''))
          print(type)
          stopifnot(paste("Color type is not recognized. Try something like \"gray.colors\"") == 999)
     }
     if (reverse) { col <- rev(col) }

     assert.agw(length(col) == totalColors)
     
     return(col)
}




## =================================================================
## pdf.for.heatmap.agw:
## * This is used INSTEAD of a call to pdf
## * Automatically calculates the proper PDF dimentions for a heatmap.agw plot.
## * This way the heatmaps can be scaled to have cells that are semi-approximately the right size.
##   Otherwise, heatmaps with lots of rows tend to have really cramped rows, and heatmaps with few rows
##   have gigantic super-tall cells.
## * You can pass in either the matrix that you are about to plot, OR the numRows you will plot.
## =================================================================
pdf.for.heatmap.agw <- function(file=file, mat=NULL, numRows=NULL, width="should not be specified by the user!!", height="should not be specified by the user!!", ...) {
     if (!is.null( mat )) {
          amount <- nrow( mat )
          assert.agw(is.null(numRows), paste("You cannot specify both a matrix AND a number of rows!! Input one OR the other! numRows was specified to be: ", numRows, ".", sep=''))
     } else if (!is.null(numRows)) {
          amount <- numRows
     }
     assert.agw(missing(width), "Hey! You should not specify width in pdf.for.heatmap.agw! It is always set to the same value.")
     assert.agw(missing(height), "Hey! You cannot manually specify height in pdf.for.heatmap.agw. It is auto-calculated!")
     assert.agw(is.numeric(amount), "Uh oh! Pdf.for.heatmap.agw is broken!")
     ALEX_HEATMAP_WIDTH_INCHES <- 15
     MAX_HEIGHT <- 80
     MIN_HEIGHT <- 12
     pdf(file=file,
         , width=ALEX_HEATMAP_WIDTH_INCHES
         , height=min(MAX_HEIGHT, max(6+amount*0.25, 12))
         , ...)
}

## =================================================================
## Alex's Heatmap-and-histogram plot
## This plotting function uses the "layout" split-up-the-canvas function to split the canvas in 3.
## For a PDF, this heatmap.agw must be at least 12 inches tall for the heatmap AND the histogram to both fit.
## It can be 8 or more inches wide and look OK.
## Use "pdf.for.heatmap.agw" to generate the heatmap, or you'll get out-of-bounds errors probably!
## =================================================================
heatmap.agw <- function(mmm, breaks=12, labRow=NULL, labCol=NULL, colorVec=NULL, colorStyle="heat"
                        , main="Heatmap", title="", cexRow=NULL, cexCol=1.5, maxNumLabels=1000
                        , col.names=NULL, row.names=NULL, cluster.rows=FALSE) {
     ## mmm: a matrix to plot
     ## Breaks: the number of histogram breaks, used for the color scheme
     ## labRow / labCol are row/column labels. You can probably override them... not sure. It may also default to the names from mmm?
     ## maxNumLabels: do not print labels if there are more than this many labels ***with actual non-blank content***
     ## ** colorStyle is recommended for easily changing the color scheme. It accepts anything that "colors.agw" takes:
     ##        colorStyle="greenwhitered" or "greenblackred" or "blueblackyellow" or "blueblackyellow2"
     ##                   or "gray" or "brown" or "sepia" or "heat"

     ## command for testing:
     ## NN = 40; source("~/TimeForScience/Lab_Code/R/AGW/agwUtil.R") ; m=matrix(rnorm(NN*10),ncol=NN); m[1:10,1] <- 50; m[1:5,2] <- -10; pdf.for.heatmap.agw("zog.pdf", m); heatmap.agw(m, breaks=seq(-6,6,length.out=12), colorStyle="blueBlackYellow2") ; dev.off()
     
     ## colorVec: the SPECIFIC list of colors to pass in. Must be equal in length to (breaks - 1)
     ## colorStyle: OR you can specify the colors as a style. This is one of the strings accepted by "colors.agw"--for example, "sepia" or "gray" or "blueblackyellow"
     min.raw <- min(mmm, na.rm=TRUE)
     max.raw <- max(mmm, na.rm=TRUE)
     mean.raw <- mean(mmm, na.rm=TRUE)
     median.raw <- median(mmm, na.rm=TRUE)
     
     print.agw("heatmap.agw: Now generating a \"heatmap.agw\" figure. If you get a \"figure region too large\" error,")
     print.agw("             that means your PDF/PNG wasn't big enough to fit the heatmap. This can be solved by using")
     print.agw("             a bigger pdf width or by using \"pdf.for.heatmap.agw\" to auto-compute the bounds.")
     if (is.logical(cluster.rows) && cluster.rows) {
          ## If cluster.rows is true, then we will CLUSTER the rows, kind of like how regular built-in "heatmap" does it.
          theDist <- stats::dist(mmm, method="euclidean")
          maxNonNullDist <- max(theDist, na.rm=T)
          numNA <- sum(is.na(theDist))
          if (numNA > 0) {
               warnMsg <- "Note, in heatmap.agw, there are 'NA' values in the distance matrix! This is sub-optimal, as stats::hclust can't actually handle NA values. Thus, we have an inelegant workaround: the NA values for distance are being replaced by the maximum non-NA value in the heatmap, for purposes of clustering. This usually results in an acceptable result."
               warning(warnMsg) ; log.red.agw(warnMsg)
               theDist[is.na(theDist)] <- maxNonNullDist ## replace any NA values with the maximum non-NA distance
          }
          hcc <- stats::hclust(theDist, method="complete") # Note: hclust can't handle NA values!
          dcc <- as.dendrogram(hcc)
          reordering <- order.dendrogram(dcc)
          mmm <- mmm[reordering, , drop=F] ## REORDER THE INPUT MATRIX BASED ON THE CLUSTERING
          assert.agw(is.null(labRow) && is.null(row.names), "Uh oh! You cannot specify that you want the matrix to be re-clustered AND ALSO specify row names. This is because once you recluster, the row names will not be what you probably expect! i.e., the row names move around!")
     }
     
     if (is.vector(mmm)) {
          mmm <- as.matrix(mmm) ## Turn a vector into a one-column-or-something matrix, so we can assume it's a matrix below.
     }
     stopifnot(is.matrix(mmm)); stopifnot(nrow(mmm) >= 1); stopifnot(ncol(mmm) >= 1)
     
     # Generates a three-row multi-part figure.
     # Top part: the caption (title)
     # Middle part: the histogram and distribution key (probably should be made optional, but it's required for now)
     # Bottom part: the heatmap
     
     if (!missing(col.names) && !is.null(col.names)) { labCol <- col.names } ## "col.names" is just an alias for "labCol"
     if (!missing(row.names) && !is.null(col.names)) { labRow <- row.names } ## "row.names" is just an alias for "labRow"
     if (is.null(labRow)) { labRow = rownames(mmm) }
     if (is.null(labCol)) { labCol = colnames(mmm) }
     
     scale01 <- function(x, low = min(x), high = max(x)) { return((x - low)/(high - low)) } ## Scale whatever the previous range was, now from 0 to 1. So like, -49 to 738 would be rescaled 0 to 1 where -49 would now be zero, and 738 would now be 1.

     if (missing(breaks) || is.null(breaks) || (length(breaks) == 1)) {
          ## If breaks was not specified, OR it was a length-one scalar, then make it into a sequence
          if (!is.null(breaks) && (length(breaks) == 1) && breaks < 3) {
               warning(paste("If you specify the number of breaks as a SINGLE NUMBER (not a sequence of numbers / vector), then that number must be >= 3. Your specified 'breaks' number was only ", breaks, ".", sep=''))
               stopifnot(breaks >= 3)
          }
          breaks.vec <- seq(min.raw, max.raw, length.out=breaks)
     } else {
          stopifnot(is.vector(breaks)) # "If 'breaks' is specified here, it must be a vector or single number."
          breaks.vec = breaks
     }
     stopifnot(length(breaks.vec) >= 3) # breaks.vec needs to be at least 3 elements long by this point!

     numValuesOutOfBoundsBelow <- sum(mmm < min(breaks.vec), na.rm=T)
     numValuesOutOfBoundsAbove <- sum(mmm > max(breaks.vec), na.rm=T)
     
     layout(matrix(c(1,2,3), byrow=T, ncol=1), heights=c(lcm(6),lcm(12),1) ) # <-- we plot THREE things, stacked vertically. This splits up the canvas into three sub-plots.
     # ==========================================
     # ======================================

     outOfBoundsText = ''
     if (numValuesOutOfBoundsAbove > 0) {
          outOfBoundsText <- paste(outOfBoundsText, numValuesOutOfBoundsAbove, " values greater than ", max(breaks.vec), " are shown in the maximum bin in the histogram. ", sep='')
     }
     if (numValuesOutOfBoundsBelow > 0) {
          outOfBoundsText <- paste(outOfBoundsText, numValuesOutOfBoundsBelow, " values less than ", min(breaks.vec), " are shown in the minimum bin in the histogram. ", sep='')
     }
     if (numValuesOutOfBoundsAbove > 0 || numValuesOutOfBoundsBelow > 0) {
          outOfBoundsText <- paste(outOfBoundsText, "\n", sep='') ## Add a newline
     }
     ## This is the FIRST of three "layout" sub-items (which are stacked vertically). It's a text description of what is being plotted.
     textDescriptionPlotAGW(paste("Heatmap with ", length(mmm)
                                  , " values, in ", nrow(mmm)
                                  , " rows and "
                                  , ncol(mmm), " columns."
                                  , " Mean = ", format(mean.raw, digits=3, nsmall=1)
                                  , ", median = ", format(median.raw, digits=3, nsmall=1)
                                  , ". "
                                  , outOfBoundsText
                                  , title
                                  , sep=''), wraplen=150, leftMargin=0.05, topMargin=0.05, cex=1.4)
     # ======================================
     # ==========================================
     # ==========================================
     # ======================================

     ## This is the SECOND of three layout things. It's a histogram
     par(mar=c("bottom"=3, "left"=9.5, "top"=5, "right"=15.5))
     
     if (missing(colorVec) || is.null(colorVec) || (length(colorVec) == 0)) {
          assert.agw(length(colorStyle) == 1, "colorStyle needs to be a string like blueblackyellow. The user can pass in a string here to automagically pick the color scheme. Or they can specify it manually with \"colorVec\".")
          ## colorStyle is a CHARACTER vector. Options include "greenwhitered" "blueblackyellow" and "blueblackyellow2" and "gray" and "sepia" . "heat" is also popular.
          colorVec <- colors.agw(n=(length(breaks.vec)-1), type=colorStyle)
     }
     
     ## Draw the background for the "legend" histogram
     COLOR_ZLIM <- c(min(breaks.vec), max(breaks.vec)) ## zlim has something to do with the maximum/minimum colors or something.
     ZERO_TO_ONE_SCALE_VEC <- c(0,1,0,1) # Scales the histogram (and image!) to fit into a 0-to-1 x and y axis scale. So the left is 0, and the right is 1.0.
     ZERO_TO_ONE_PLUS_EXTRA_ON_Y_AXIS_SCALE_VEC <- c(0,1,-0.05,1.05) # The scale should be SLIGHTLY different from the ZERO_TO_ONE_SCALE_VEC--it needs to go slightly lower and slightly higher, so as to not clip the values off the bottom & top of the histogram. Thus, instead of 0 and 1, we use -0.05 and 1.05
     ## ============== DRAW LEGEND BACKGROUND GRADIENT =================
     par(usr=ZERO_TO_ONE_SCALE_VEC) # Scales the histogram (and image!) to fit into a 0-to-1 x and y axis scale. So the left is 0, and the right is 1.0.
     legendBackground.mat <- matrix(seq(min(breaks.vec), max(breaks.vec), length.out=length(colorVec)), ncol=1)
     image(z=legendBackground.mat, col=colorVec, breaks=breaks.vec, zlim=COLOR_ZLIM, xaxt="n", yaxt="n") ## This is the histogram / distribution background that goes at the top of the heatmap.

     ## ============== DRAW THE 'HERE IS THE MEDIAN' VERTICAL LINE =================
     MEDIAN_LINE_WIDTH    <- 3
     abline(v=scale01(median.raw, min(breaks.vec), max(breaks.vec)), lwd=2*MEDIAN_LINE_WIDTH, lty="solid", col="gray")
     abline(v=scale01(median.raw, min(breaks.vec), max(breaks.vec)), lwd=MEDIAN_LINE_WIDTH, lty="dashed", col="black")
     
     ## ============== DRAW THE WIGGLY HISTOGRAM LINE =================
     par(usr=ZERO_TO_ONE_PLUS_EXTRA_ON_Y_AXIS_SCALE_VEC) # The scale should be SLIGHTLY different from the ZERO_TO_ONE_SCALE_VEC--it needs to go slightly lower and slightly higher, so as to not clip the values off the bottom & top of the histogram. Thus, instead of 0 and 1, we use -0.05 and 1.05
     HISTOGRAM_LINE_WIDTH <- 6

     mmmClippedToBounds <- mmm
     mmmClippedToBounds[ mmmClippedToBounds < min(breaks.vec) ] <- min(breaks.vec)
     mmmClippedToBounds[ mmmClippedToBounds > max(breaks.vec) ] <- max(breaks.vec)
     
     histValues <- hist(mmmClippedToBounds, plot=FALSE, breaks=breaks.vec) # XLIM is not used except when plotting! xlim=c(min(breaks.vec)*1.05, max(breaks.vec)*1.05))
     hx <- scale01(breaks.vec, min(breaks.vec), max(breaks.vec))
     hy <- c(histValues$counts, histValues$counts[length(histValues$counts)])
     scaledHy <- hy/max(hy)
     lines(hx, scaledHy, lwd=2*HISTOGRAM_LINE_WIDTH, type="s", col="white") ## Draw the actual histogram as a squiggly line
     lines(hx, scaledHy, lwd=HISTOGRAM_LINE_WIDTH, type="s", col="black")
     ## ============== DRAW LEGEND AXES LABELS FOR THE HISTOGRAM / LEGEND =================
     labelsPretty <- pretty(breaks.vec)
     axis(1, at=scale01(as.numeric(labelsPretty), min(breaks.vec), max(breaks.vec)), labels=round(labelsPretty, 2)) ## X AXIS
     axis(2, at=pretty(hy)/max(hy), pretty(hy))
     ## ============== DRAW BOX AROUND THE HISTOGRAM / LEGEND =================
     box(lwd=1)
     ## ============== DRAW TITLE =================
     par(cex.main=1.7) # No idea why we scale it to 1.7, but I guess that is a decent size
     title(paste("Color Key and Histogram of the Distribution of Heatmap Values\nDashed line at median value (", format(median.raw, digits=3, nsmall=1), ")", sep=''))
     ## ============== DRAW INFORMATIONAL TEXT AT VERY TOP OF GRAPH =================
     mtext(side=2, paste("Count (out of ", length(mmmClippedToBounds) ," total)", sep=''), line=3)
     # ==========================================
     # ======================================

     ## This is the THIRD thing being plotted, and is the real heatmap. If you want a heatmap that doesn't use layout/mfrow/etc., you can just copy out this code here.
     ## This is the "meat" of the heatmap generation.
     
     par(mar=c("bottom"=25, "left"=8, "top"=8, "right"=20), cex.main=2.2, cex=0.5)
     mainHeatmap.mat <- t(mmmClippedToBounds)[, nrow(mmmClippedToBounds):1, drop=F]  ## <-- transpose AND flip, to rotate the correct way
     mainWithDimensions.string <- paste(main, "\n(", nrow(mmmClippedToBounds), " rows by ", ncol(mmmClippedToBounds), " columns)", sep='')
     image(mainHeatmap.mat, breaks=breaks.vec, axes=F, main=mainWithDimensions.string, col=colorVec, zlim=COLOR_ZLIM)
     # Notice that ‘image’ interprets the matrix as a table of
     # ‘f(x[i], y[j])’ values, so that the x axis corresponds to row
     # number and the y axis to column number, with column 1 at the
     # bottom, i.e. a 90 degree counter-clockwise rotation of the
     # conventional printed layout of a matrix.

     numNonBlankRows <- 0
     if (!is.null(labRow)) {
          numNonBlankRows <- sum(!is.na(labRow) & (nchar(labRow) > 0), na.rm=T) ## Count the number of NON-BLANK rows only!
     }

     if (missing(cexCol) || is.null(cexCol)) {
          cexCol <- min(2.0, 40.0/ncol(mmm))
     }
     if (is.null(cexRow)) {
          rowNumberToUseForSizeCalculation <- nrow(mmm) ## Make the labels tiny enough to individually specify a single row
          cexRow = 1.0
          rowsAtWhichCexIsMin <- 800 ## The number of rows at which the text is the tiniest (or more rows than this)
          rowsAtWhichCexIsMax <- 100 ## The number of rows at which the text is the biggest (or fewer rows than this)
          stopifnot(rowsAtWhichCexIsMin > rowsAtWhichCexIsMax)
          maxRowCex <- 1.00
          minRowCex <- 0.10

          if (rowNumberToUseForSizeCalculation >= rowsAtWhichCexIsMin) {
               cexRow <- minRowCex
          } else if (rowNumberToUseForSizeCalculation <= rowsAtWhichCexIsMax) {
               cexRow <- maxRowCex
          } else {
               fractionTowardMinimum <- (rowNumberToUseForSizeCalculation - rowsAtWhichCexIsMax)/(rowsAtWhichCexIsMin - rowsAtWhichCexIsMax)
               cexRow <- 0.10 + (maxRowCex-minRowCex)*(1 - (fractionTowardMinimum**0.5))
          }
     }

     if (!is.null(labRow) && numNonBlankRows <= maxNumLabels) {
          ## Right side axis--got to FLIP them around, because of the way we printed the image
          tickLoc.vec <- (0:(nrow(mmm)-1))/(nrow(mmm)-1)
          backwardsLabels.vec <- rev(labRow) ## Gotta reverse the order here!
          ## You can verify that this is correct by examining the values in "heatByCoef" and verifying that
          ## they match up here, too.
          axis(4, at=tickLoc.vec, labels=backwardsLabels.vec, tick=F, las=2, cex.axis=cexRow) # 4 = right axis, usually with gene names
     }

     if (ncol(mmm) == 1) {
          bottomAxisLoc.vec <- c(0.5) ## special case for ncol == 1, so we don't divide by zero!
     } else {
          bottomAxisLoc.vec <- seq(from=0, to=(ncol(mmm)-1)) / (ncol(mmm)-1)
     }
     
     axis(1, at=bottomAxisLoc.vec, labels=labCol, tick=F, las=2, cex.axis=cexCol) # 1 = bottom axis, usually with array names
     box(lwd=1)

     print.agw("heatmap.agw: Looks like the figure region was an acceptable size to fit the heatmap.")
     # ======================================
     # ==========================================
}








## =============================================================================
## plot.multi.agw: a kind of hack-ish way to plot multiple types of graph from
## the same call. Like if you wanted both a high-res PDF and also a PNG of the same data.
## This is probably not a particularly useful function in most/all cases.

## Example call:

##         plot.multi.agw(type.func=c(
##                        call("pdf", paste(heatmapPath, ".pdf", sep=''), width=15, height=20)
##                        , call("png", paste(heatmapPath, ".png", sep=''), pointsize=18, res=72, width=1800, height=2400))
##                        , plot.func=call("heatmapActual.function")
##                        )
## =============================================================================
plot.multi.agw <- function(type.func, plot.func) {
     ## The arguments to this must be functions generated with "call"
     ## See the examples in the assertions below.
     
     ## type.func needs to be something like "png(filename)" or "pdf(filename)"
     ## plot.func the actual how-to-make-a plot function (example: plot / hist)
     
     for (thePlotMethod in type.func) {
          assert.agw(is.language(thePlotMethod), "The type.func argument must be a VECTOR of types of R functions to output plot data. Example: type.func=c(call(\"png\", width=800), call(\"pdf\", width=10, height=7))")
          assert.agw(is.language(plot.func), "The plot.func argument is the *single* function to call to actually generate the plot. It must be one function, generated with \"call.\" Example: plot.func=call(\"plot\", 1:10, col=\"blue\");")
          print.cyan.agw("plot.multi.agw is running...")
          eval.parent(thePlotMethod) # <-- evaluate in the CALLING environment
          eval.parent(plot.func)     # <-- evaluate in the CALLING environment
          dev.off()
     }
     
}




### ===============================================================================
## Input: a matrix (mm)
## Now, this function goes through each column of the matrix and plots the data from that column as a line.
## If you want to do it by rows instead, transpose the input matrix.
## col can be a VECTOR of colors. So you can pass in something like rainbow( ncol(mm) ) if you want.
## What you end up with is a nice plot with a bunch of lines, like you might see on a sales chart of
## various car brands over time, with several wiggly lines.
## Note that you don't specify the X-axis: it is auto-set to be one integer per Y data point.

## This appears to be basically a reimplementation of "matplot," so you should probably
## use that instead.

### ===============================================================================
## plotLinesFromMatrixColumnsAgw <- function(mm=NULL, ..., ylim=NULL, col="black", lwd=1, lty=1) {

##      ## This is probably a duplication of the built-in function "matplot"!

##      warning("I don't think plotLinesFromMatrixColumnsAgw is probably a good function---you should probably just use the built-in function <matplot> instead!!")
##      if (!is.matrix(mm) && is.vector(mm) && (length(mm) > 0)) {
##           warning("plotLinesFromMatrixColumnsAgw: the input matrix was actually a one-dimensional vector! Trying to work around this...")
##           mm <- matrix( mm , byrow=FALSE, ncol=1)
##      }
##      #assert.agw(!is.null(mm) && is.matrix(mm), "Problem! Input to plotLinesFromMatrixColumnsAgw must be a non-null matrix!")
##      if (length(col) < ncol(mm)) { col <- c(col, rep("#00000088", times=ncol(mm))) } ## default fallback color is semi-transparent black
##      plot(  x=(1:nrow(mm))
##           , y=mm[,1, drop=TRUE]
##           , type='l'
##           , ylim=ylim
##           , ...
##           , lwd=lwd, lty=lty, col=col[1])
##      for (i in 2:ncol(mm)) {
##           lines(  x=(1:nrow(mm))
##                 , y=(mm[,i, drop=TRUE]), lwd=lwd, lty=lty, col=col[i])
##      }
## }
## ### ===============================================================================
## ### ===============================================================================










### ===============================================================================
### Used in the scatterplots that compare values across two arrays
### ===============================================================================
panel.correlation.local <- function(x, y, digits=2, prefix="", cex.cor=4.0, boxWidth=2) {
    usr <- par("usr"); on.exit(par(usr)) ## restore settings on finishing the plot
    SQUARE_MAX <- 1.0 ; SQUARE_MIN <- 0.0
    MIDDLE_OF_SQUARE <- (SQUARE_MAX + SQUARE_MIN)/2.0
    par(usr = c(SQUARE_MIN, SQUARE_MAX, SQUARE_MIN, SQUARE_MAX))
    #‘usr’ A vector of the form ‘c(x1, x2, y1, y2)’ giving the extremes of the user coordinates of the plotting region.
     
    origR  <- cor(x, y, use="pairwise.complete.obs", method="pearson")
    absR   <- abs(origR)
    MIN_CEX_FAC <- 0.4

	theCex <- cex.cor * max(MIN_CEX_FAC, absR**2, na.rm=T) ## don't let the CEX get any smaller than the MIN_CEX

    backgroundColor <- hsv(s = 0.0,  v = 1.0 - 5*(1.0 - max(0.9, absR, na.rm=T))) ## uhh... 1.0 = white, and lower scores = a dark-ish gray... kind of hack-ish. Basically it gets as dark as it's going to get around R=0.7 or thereabouts.
    rect(xleft=-100, ybottom=-100, xright=100, ytop=100, col=backgroundColor)
    text(MIDDLE_OF_SQUARE, MIDDLE_OF_SQUARE
         , format(origR, digits=2, nsmall=2, scientific=FALSE)
         , cex=theCex ## <-- The closer R is to zero, the smaller the text size
         )
    box(lwd=boxWidth)
}

### ===============================================================================
## Makes a big multi-plot matrix of correlations between arrays / values
## Uses the data in "dataMatrix" to figure out what to plot.
## if you for some reason actually want NO labels, set labelVec to "" in the arguments
### ===============================================================================
pairsCorMatrixPlotAGW <- function(filePath, dataMatrix, labelVec=NULL, keys, main="log2(Intensity).  Red points = within-group comparison. Values in lower left are Pearson's R of the log-transformed values.") {
     print("Note: this is an OLD function for microarrays only, you should use pairs.agw for all future code.")
     print("Note: this is an OLD function for microarrays only, you should use pairs.agw for all future code.")
     print("Note: this is an OLD function for microarrays only, you should use pairs.agw for all future code.")
     warning("Note: this is an OLD function for microarrays only, you should use pairs.agw for all future code.")
     assert.agw(nrow(keys$"table") == ncol(dataMatrix), "problem! keys were not the same length as the dataMatrix.")
     if (missing(labelVec) || is.null(labelVec)) {
          COLUMN_THAT_PROBABLY_HAS_THE_FILENAMES <- 1 ## usually the first column in the keys file...
          filenamesMostLikely.vec <- keys$table[,COLUMN_THAT_PROBABLY_HAS_THE_FILENAMES]
          labelVec <- paste(filenamesMostLikely.vec, "\n(", keys$groups[keys$group.assignments], ")", sep='')
     }
     labelVec <- gsub("[ ]+", "\n", labelVec, perl=TRUE, ignore.case=TRUE) ## labelVec with a newline added for any spaces
     labelVec <- gsub("\\.(CEL|gpr)", " ", labelVec, perl=TRUE, ignore.case=TRUE )  ## put a newline before the .CEL or .gpr extension at the end of the filename
     labelVec <- gsub("[ _]+vs[.]?[ _]+", "\nvs.\n", labelVec, perl=TRUE, ignore.case=TRUE) ## break up "vs" so it spans multiple lines
     labelVec <- gsub("\n\n+", "\n", labelVec, perl=TRUE, ignore.case=TRUE) ## One newline in a row at most!
     labelVecNoNewlines <- gsub("\n", "", labelVec, perl=TRUE, ignore.case=TRUE)
     print.agw("Drawing a \"pairs\" plot matrix with many sub-plots to file <", filePath, ">... (note, this can take a couple of minutes)\n")
     ## cumul/whichBin: Shows which "bins" each number is in. So like, experiment #2 might be in experimental group #1, if there were 2 replicates of that experiment.

     exprs.in.each.group <- table(keys$group.assignments)
     cumul <- c(0, cumsum(exprs.in.each.group)) ## <-- 0 is the leftmost boundary!
     whichBin <- function(value, binMarkers) { ## right sides! The first value is NOT a valid one
          ## Tells you which "bin" an array is in. The bin is the experimental group.
          ## So if you have 9 arrays, in 3 groups (control, control, control, exp1, exp1, exp1, exp2, exp2, exp2),
          ## then the bins are (1, 1, 1, 2, 2, 2, 3, 3, 3)
          ## binMarkers must be sorted in ascending order!
          for (i in 1:(length(binMarkers)-1)) {
               left  <- binMarkers[i]   ;  right <- binMarkers[i+1]
               if (value > left && value <= right) { return(i); }
          }
          return(-1); ## out of bounds
     }

     nGroups    <- length(keys$"groups")
     pointAlpha <- 0.75  ## 0.75 = 75% opaque.
     regularPointColor <- hsv(h=1, s=1, v=0, alpha=pointAlpha) ## The foreground color for points in each "bin"
     backgroundAlpha   <- 1.0
     backColors        <- rainbow(n=nGroups*(nGroups-1), s = 0.25, v = 1.0, alpha=backgroundAlpha) ## The background colors for each "bin"
     WITHIN_GROUP_POINT_COLOR <- hsv(h=1, s=1, v=1, alpha=pointAlpha)
     WITHIN_GROUP_BACKGROUND_COLOR   <- "white" ## What color are the cells that have the array filenames?
     MINI_PLOT_BORDER_THICKNESS <- 2.0   ## How thick the border around each scatterplot is

     if (ncol(dataMatrix) > 30) {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 100
     } else if (ncol(dataMatrix) > 60) {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 75
     } else {
          SIZE_IN_PIXELS_PER_SCATTERPLOT <- 200 ## Size of each mini-scatterplot
     }
     MIN_ALLOWED_SIZE_OF_PLOT_IN_PIXELS     <- 1000 ## Minimum size of this plot
     TOTAL_PNG_SQUARE_SIZE <- max((ncol(dataMatrix) * SIZE_IN_PIXELS_PER_SCATTERPLOT), MIN_ALLOWED_SIZE_OF_PLOT_IN_PIXELS)
     # ======================================

     png(filePath, height=TOTAL_PNG_SQUARE_SIZE, width=TOTAL_PNG_SQUARE_SIZE, pointsize=18, res=72)
     par(mar=c("bottom"=2, "left"=2, "right"=2, "top"=8))
     par(cex.axis=1.0, cex.lab=1.2, cex.main=1.0)     
     dimSize <- ncol(dataMatrix)
     whichSubplot <- list("x" = 2, "y" = 1)  ## <-- have to MANUALLY keep track of which sub-plot we are plotting in the  big "upper.panel" function using these vars and the <<- assignment operator
     corrMatrixOutputFilePath <- gsub(".png", ".pearson.corr.matrix.out.txt", filePath)
     print.green.agw("Now writing output correlations to <", corrMatrixOutputFilePath,">...")
     vvv <- matrix(nrow=length(labelVecNoNewlines), ncol=length(labelVecNoNewlines), dimnames=list(labelVecNoNewlines, labelVecNoNewlines))
     for (xx in seq_along(labelVecNoNewlines)) {
          for (yy in seq_along(labelVecNoNewlines)) {
               vvv[xx, yy] <- cor(dataMatrix[,xx], dataMatrix[,yy], use="complete", method="pearson")
          }
     }
     NUM_DECIMAL_PLACES_FOR_CORR_OUTPUT <- 4
     write.table(signif(vvv, NUM_DECIMAL_PLACES_FOR_CORR_OUTPUT), corrMatrixOutputFilePath, sep="\t", col.names=NA, row.names=TRUE, quote=T)
     print.green.agw("[Done] writing output correlations to <", corrMatrixOutputFilePath,">.")
     
     ## "graphics::pairs" plot is one of those big matrix plots that is ususally used to display many scatterplots at once. In this case, it's scatterplots (top right -- upper.panel) and correlation values (bottom left -- lower.panel)
     graphics::pairs(dataMatrix
                     , main=main, labels=labelVec
                     , cex.labels=1.5 ## <-- this determines the labels on the diagonal
                     , cex.main=1.0 ## <-- this is related to the size of the graph in a non-intuitive way
                     , lower.panel=function(x,y) {
                          #rect(-10,-10,10,10, col="gray", border=NA) ;
                          panel.correlation.local(x,y, boxWidth=MINI_PLOT_BORDER_THICKNESS) ;
                     }
                     , upper.panel=function(x,y) {
                          binX <- whichBin(whichSubplot$x, cumul)
                          binY <- whichBin(whichSubplot$y, cumul)
                          binQuasiRandomIdentifier <- (((binY-1) * nGroups * 13) + binX*17)  ## has to be somehow based on the bins, so that everything in the same bin will be the same color. The exact formula is irrelevant, though, as long as it "seems" random.
                          if (binX == binY) {
                               thePointColor <- WITHIN_GROUP_POINT_COLOR
                               theBackColor <- WITHIN_GROUP_BACKGROUND_COLOR
                          } else {
                               thePointColor <- regularPointColor
                               theBackColor <- backColors[1 + (binQuasiRandomIdentifier %% length(backColors))]
                          }
                          cat(paste("Now calculating the box at X,Y (", whichSubplot$x, ", ", whichSubplot$y, "), which is in bin (", binX, ", ", binY, "), numbered ", binQuasiRandomIdentifier, " and colored <", theBackColor, ">... (", length(x), " points)", sep=''))
                          rect(xleft=min(x)-10, ybottom=min(y)-10, xright=max(x)+10, ytop=max(y)+10, col=theBackColor, border=NA)
                          points(x,y, pch='.', cex=4, col=thePointColor) ## <-- pch='.' is 10 times faster than any other plot character. Use cex=### to adjust the size of the dot.
                          box(lwd=MINI_PLOT_BORDER_THICKNESS)
                          cat(" [Done]\n")
                          whichSubplot$x <<- (whichSubplot$x+1); ## go to the next column... note that this is a GLOBAL assignment! That is important.
                          if (whichSubplot$x > dimSize) {
                               whichSubplot$x <<- (whichSubplot$y+2) ## <-- because upper.panel is a right triangle, we start from a point farther to the right on the plot each time. (This really is supposed to be (xVal <<- yVal+2)!
                               whichSubplot$y <<- (whichSubplot$y+1) ## <-- this only works because the upper.panel is a *right triangle*, so it just happens that the X axis index where we start a new row is the same as the Y axis index. Y axis index is also the same as the X##(whichSubplot$y) #startX;
                          }
                     })
     dev.off()
     print.agw("[Done] with pairs plot.")
}

### ===============================================================================
### pairs.agw
### ===============================================================================

pairs.agw <- function(data, groups.vec=NULL, main="Pairs plot. Red points = within-group comparison. Values in lower left are Pearson's R.", bw=FALSE) { # plot pairs
     # "data": a (numeric) data matrix.
     # "groups.vec": a vector of the experimental groups for each COLUMN in the data matrix: example: groups.vec=c("WT","WT","DrugX","DrugY")
	# bw (default: FALSE). Set to true to make the output black and white only. Useful for some publication figures.
	
     labelVec <- colnames(data)
     if (is.null(labelVec)) { labelVec <- paste("Column_", 1:ncol(data), sep='') }

     if (is.null(groups.vec)) {
          groups.vec <- 1:ncol(data)  # Default: puts every item into its own group, but DON'T bother listing these (mostly useless) numeric groups in the labelVec.
     } else {
          stopifnot(is.vector(groups.vec)); stopifnot(ncol(data) == length(groups.vec)) # Every array (column) must be in a group, if groups were manually specified.
          labelVec <- paste(labelVec, "\n(", groups.vec, ")", sep='') # <-- only add groups to the name if they were in fact specified in the function call
     }
     binAssignments.vec <- as.numeric(as.factor(groups.vec))
     
     themin <- min(data,na.rm=T); themax=max(data,na.rm=T); mmdist <- themax-themin
     dataWithTwoFakeRowsForPlotting <- rbind(data, "MIN_ROW_AT_END"=rep(themin, times=ncol(data)), "MAX_ROW_AT_END"=rep(themax, times=ncol(data)) )
     
     nGroups    <- nlevels(as.factor(groups.vec))
     pointAlpha <- 0.50  ## 0.50 = 50% opaque.

     if (bw) {
     	# Black and white only
	regularPointColor <- hsv(h=1, s=0, v=0, alpha=pointAlpha) # black
	backgroundAlpha <- 0.0
        WITHIN_GROUP_BACK_COLOR   <- "#FFFFFF00" # transparent
        WITHIN_GROUP_POINT_COLOR <- hsv(h=1, s=0, v=1, alpha=pointAlpha) # black
     } else {
     	# Colors please
	regularPointColor <- hsv(h=1, s=1, v=0, alpha=pointAlpha) ## The foreground color for points in each "bin"
	backgroundAlpha   <- 1.0
        WITHIN_GROUP_BACK_COLOR   <- "white"
	WITHIN_GROUP_POINT_COLOR <- hsv(h=1, s=1, v=1, alpha=pointAlpha) # red
     }

     allBackColors     <- rainbow(n=nGroups*(nGroups-1), s = 0.25, v = 1.0, alpha=backgroundAlpha) ## The background colors for each "bin"
     PLOT_BORDER_THICKNESS <- 1.0   ## Thickness of the border around each scatterplot
     par(mar=c("bottom"=2, "left"=2, "right"=2, "top"=8))
     par(cex.axis=1.0, cex.lab=1.2, cex.main=1.0)     
     dimSize <- ncol(dataWithTwoFakeRowsForPlotting)
     whichSubplot <- list("x"=2,"y"=1)  ## <-- have to MANUALLY keep track of which sub-plot we are plotting in the "upper.panel" function using these vars and the <<- global assignment operator

     panel.correlation.agw <- function(x, y, digits=2, cex.cor=4.0, boxWidth=2) {
          usr <- par("usr"); on.exit(par(usr)) ## restore settings on finishing the plot
          par(usr=c(0,1,0,1)) # <- gives th extremes of the user coordinates of the plotting region.
          # Note: correlation is NOT computed on the last two elements: (hence the x-2, y-2 below). This is because those are FAKE maximum / minimum points that are required in order to make graphics::pairs plot everything on the same scale.
          origR   <- cor(x[1:(length(x)-2)], y[1:(length(y)-2)], use="pairwise.complete.obs", method="pearson")
          absR   <- abs(origR)
          MIN_CEX_FAC <- 0.4
          corCex <- cex.cor*max(MIN_CEX_FAC, absR**2, na.rm=T) ## don't let the CEX get any smaller than the MIN_CEX. Correlations closer to zero have smaller text size.
	  if (bw) {
	  	corColor <- "black"
		backgroundColor <- "#FFFFFF00" # white, fully transparent
	  } else {
	  	  if (origR < 0) { corColor = "#990000" } else { corColor = "black" } # Correlation color (negative is dark RED)
	          backgroundColor <- hsv(s=0.0,  v=1.0 - 5*(1.0 - max(0.9, absR, na.rm=T))) ## v=1.0 is white, and lower scores = a dark-ish gray. Basically it gets as dark as it's going to get around R=0.7 or thereabouts.
          }
          rect(xleft=-1, ybottom=-1, xright=2, ytop=2, col=backgroundColor)
          text(0.5, 0.5, paste("r=", format(origR, digits=digits, nsmall=digits, scientific=FALSE), sep='')
               , cex=corCex, col=corColor) ## <-- The closer the correlation value is to zero, the smaller the text size
          box(lwd=boxWidth)
     }

     par(pty='s') # Square plot
     graphics::pairs(dataWithTwoFakeRowsForPlotting
                     , main=main
                     , labels=labelVec
                     , cex.labels=1.5 ## <-- this determines the labels on the diagonal
                     , cex.main=1.0 ## <-- this is related to the size of the graph in a non-intuitive way
                     , lower.panel=function(x,y) {
                          #rect(-10,-10,10,10, col="gray", border=NA) ;
                          panel.correlation.agw(x,y, boxWidth=PLOT_BORDER_THICKNESS) ;
                     }
                     , upper.panel=function(x,y) {
                          binX <- binAssignments.vec[whichSubplot$x]
                          binY <- binAssignments.vec[whichSubplot$y]
                          binQuasiRandom <- (((binY-1) * nGroups * 13) + binX*17)  ## hash it to get a "random" color
                          pointColor <- ifelse((binX == binY), WITHIN_GROUP_POINT_COLOR, regularPointColor)
                          backColor  <- ifelse((binX == binY), WITHIN_GROUP_BACK_COLOR , allBackColors[1 + (binQuasiRandom %% length(allBackColors))] )
                          cat(paste("pairs.agw: calculating the box at (", whichSubplot$x, ", ", whichSubplot$y, ") in bin (", binX, ", ", binY, ") with ", length(x), " points.\n", sep=''))
                          rect(xleft=themin-mmdist, ybottom=themin-mmdist, xright=themax+mmdist, ytop=themax+mmdist, col=backColor, border=NA)
                          points(x[1:(length(x)-2)],y[1:(length(y)-2)], pch='.', cex=4, col=pointColor) ## <-- pch='.' is 10x faster than any other plot character
                          # Note: we do NOT plot the last two elements (hence the x-2, y-2 above). This is because those are FAKE maximum / minimum points that are required in order to make graphics::pairs plot everything on the same scale.
                          abline(a=0, b=1, lty=2, lwd=1, col="black")
                          box(lwd=PLOT_BORDER_THICKNESS)
                          whichSubplot$x <<- (whichSubplot$x+1); ## next column... note that this is a GLOBAL assignment (<<-)
                          if (whichSubplot$x > dimSize) {
                               whichSubplot$x <<- (whichSubplot$y+2) ## <-- upper.panel is a right triangle, so we start from a point farther to the right on the plot each time.
                               whichSubplot$y <<- (whichSubplot$y+1) ## <-- pper.panel is a *right triangle*, so it just happens that the X axis index where we start a new row is the same as the Y axis index.
                          }
                     })
}


### ===============================================================================
### 
### ===============================================================================

agwWriteDNAStringSetNoLinebreaks <- function(dss, file) { # Writes a DNAStringSet such that each sequence is only ONE line
  stopifnot("DNAStringSet" %in% class(dss)); stopifnot(is.character(file))
  if (is.null(names(dss))) { stop("Hey, you should set the 'names(...)' of your DNAStringSet, so that we can make named fasta records!") }
  cat( paste0(">", names(dss), "\n", as.character(dss))   , file=file, sep="\n")
}


plot_of_character_frequency_by_position <- function(summarizedCountsFile) {
          # Makes a line plot showing the base frequency by base *position*.
          # Input is normally a "summarized counts" file, which is generated by count_char_freq_at_each_position.pl.
          # This is normally used for RNA-Seq, but can be used for any set of DNA (or other) sequences.

          # summarizedCountsFile is the table that is output by the script "count_char_freq_at_each_position.pl"
          # which is normally run on a fasta file.
          # Note: if you want to get every Nth line from a file, starting with line M:
          #    cat YOURFILE | sed -n 'N~Mp' > OUTPUTFILE
          dat <- read.delim(file=summarizedCountsFile, as.is=TRUE, check.names=FALSE, row.names=1, header=1)
          RSUM_INDEX <- which(rownames(dat) == "RSUM")
          ratios <- apply(dat, 2, function(ccc) { ccc/ccc[RSUM_INDEX] })
          ratios <- ratios[-RSUM_INDEX,]
          LINECOLORS <- rainbow(nrow(ratios))

          if (nrow(ratios) %in% c(4,5) && all(c("A","C","G","T") %in% rownames(dat))) {
               IS_ACGT_PLOT <- TRUE
               theYlab <- "Fraction of each base"
               LINECOLORS <- c("#FF0000BB","#000000BB","#0033FFBB","gray","#33CC00BB")
          } else {
               IS_ACGT_PLOT <- FALSE
               theYlab <- "Fraction of each character"
          }

          matplot(t(ratios), type='n'
                  , main=paste("Base frequency at each position in\n\"", summarizedCountsFile, "\"", sep='')
                  , xlab="Position in each sequence, from the FastQ file", ylab=theYlab)

          if (IS_ACGT_PLOT) {
               abline(h=0.25, col="#333333", lty="solid", lwd=2)
               abline(h=c(0.2,0.3),  col="#666666", lty="solid" , lwd=1)
          }
          matplot(t(ratios), type='b', col=LINECOLORS, lty="solid", cex=0.85, pch=rownames(ratios), lwd=2, add=T)
     }


### ===============================================================================
### 
### ===============================================================================





plotSplot <- function(x, y, main=NULL, xlab=NULL, ylab=NULL, log=TRUE, xlim=NULL, ylim=NULL, legendLocation="bottomright", legendCex=1.0) {
     ## Makes a particular kind of colored scatterplot with a bunch of blobs. I was calling it a "splot" since the colors look kind of splattered-on in abstract art fashion.
	 ## Expects the data to be passed in in NON-LOG space, but then plots the (LOG2 TRANSFORMED (value + 1.0) ) values.
     ## Adds jitter to the plot automatically.

     ## The splot is optimized for PNGs of size 1000 x 1000.
     ## You should wrap those around this call, like:
     ##   png("x.png", width=1000, height=1000)
     ##    plotSplot(...)
     ##   dev.off()

     xraw <- x
	 yraw <- y

     xIsNA.bool.vec <- is.na(xraw) # <--no value at a ALL, not even 0
     yIsNA.bool.vec <- is.na(yraw)

     EPSILON_FOR_ZERO_COUNTS <- 0.00001 ## If something is less than this, then it's close enough to zero for our purposes.
     xNoCounts.bool.vec <- (xIsNA.bool.vec | xraw < EPSILON_FOR_ZERO_COUNTS) ## NA or (close enough to) zero
     yNoCounts.bool.vec <- (yIsNA.bool.vec | yraw < EPSILON_FOR_ZERO_COUNTS) ## NA or (close enough to) zero

     xHasCounts.bool.vec <- !xNoCounts.bool.vec
     yHasCounts.bool.vec <- !yNoCounts.bool.vec

     xyNeitherHasCounts.bool.vec <- (xNoCounts.bool.vec & yNoCounts.bool.vec)
     xyExactlyOneHasCounts.bool.vec <- xor(xHasCounts.bool.vec, yHasCounts.bool.vec)
     xyBothHaveCounts.bool.vec <- (xHasCounts.bool.vec & yHasCounts.bool.vec)

     totals <- list(xTotal=length(xraw)
                    , yTotal=length(yraw)
                    , xNACount=sum(xIsNA.bool.vec)
                    , yNACount=sum(yIsNA.bool.vec)
                    , xNoneCount=sum(xNoCounts.bool.vec)
                    , yNoneCount=sum(yNoCounts.bool.vec)
                    , xHasCounts=sum(xHasCounts.bool.vec)
                    , yHasCounts=sum(yHasCounts.bool.vec)
                    , xyNoCounts=sum(xyNeitherHasCounts.bool.vec)
                    , xyOnlyOneHasCounts=sum(xyExactlyOneHasCounts.bool.vec)
                    , xyBothHaveCounts=sum(xyBothHaveCounts.bool.vec)
                    , xyHasCountInEither=sum(xyBothHaveCounts.bool.vec) + sum(xyExactlyOneHasCounts.bool.vec)
                    )

     if (log) {
          PSEUDOCOUNT_FOR_LOGSPACE <- 1
          xnew <- log2(PSEUDOCOUNT_FOR_LOGSPACE + ifelse(is.na(xraw), yes=0, no=xraw))
          ynew <- log2(PSEUDOCOUNT_FOR_LOGSPACE + ifelse(is.na(yraw), yes=0, no=yraw))
     } else {
          xnew <- ifelse(is.na(xraw), 0, xraw)
          ynew <- ifelse(is.na(yraw), 0, yraw)
     }
     
     col.vec <- rep("#FF111122", times=length(xnew))
     col.vec[xnew == 0] <- "#1100FF22"
     col.vec[ynew == 0] <- "#11FF0022"
     col.vec[xnew == 0 & ynew == 0] <- "#66666622"

     xjit=jitter(xnew, amount=0.05) ## Jittered value for display purposes
     yjit=jitter(ynew, amount=0.05) ## Jittered value for display purposes

     plot(main=main, x=xjit, y=yjit, xlab=xlab, ylab=ylab, pch=19, cex=3, col=col.vec, xlim=xlim, ylim=ylim)

     points(x=xjit, y=yjit, cex=2, pch='.', col="black")

     xZeroBecomesNA <- ifelse(xnew == 0, yes=NA, no=xnew)
     yZeroBecomesNA <- ifelse(ynew == 0, yes=NA, no=ynew)
     bestFit <- lm(yZeroBecomesNA ~ xZeroBecomesNA)
     abline(bestFit, lty=2, lwd=3, col="red")

     bothXandYareZero.bool.vec <- ((xnew == 0) & (ynew == 0))
     xOneZero <- xnew[!bothXandYareZero.bool.vec] ## Either X or Y can be zero, but not both!
     yOneZero <- ynew[!bothXandYareZero.bool.vec] ## Either X or Y can be zero, but not both!

     allCor           <- cor(x=xnew, y=ynew) ## Correlation with ALL data points (even the 0,0 ones)
     oneCanBeZeroCor  <- cor(x=xOneZero, y=yOneZero, use="all.obs")
     neitherIsZeroCor <- cor(x=xZeroBecomesNA, y=yZeroBecomesNA, use="complete.obs")

     asPercentStr <- function(numer, denom=totals$xyHasCountInEither, decimals=2) {
          return(paste(format(100*numer/denom, digits=decimals, nsmall=decimals), '%', sep=''))
     }
     legend(legendLocation, cex=legendCex, bg="#FFFFFFAA"
            , legend=c("Dashed line is the best fit line AFTER any zero values were removed"
              , paste(length(yraw), " annotated genes were examined.", sep='')
              , paste("No. genes with a non-zero value in either axis: ", totals$xyHasCountInEither, " (", asPercentStr(totals$xyHasCountInEither, denom=length(yraw)), " of all genes)", sep='')
              , paste("Below: Percentages given are relative to the total number of genes with non-zero value in either axis.", sep='')
              , paste("* > 0 value in both axes: ", totals$xyBothHaveCounts, " (", asPercentStr(totals$xyBothHaveCounts), ")", sep='')
              , paste("* > 0 value in X: ", totals$xHasCounts, " (", asPercentStr(totals$xHasCounts), ")", sep='')
              , paste("* > 0 value in Y: ", totals$yHasCounts, " (", asPercentStr(totals$yHasCounts), ")", sep='')
              , paste("* > 0 value in one method but not the other: ", totals$xyOnlyOneHasCounts, " (", asPercentStr(totals$xyOnlyOneHasCounts), ")", sep='')
              , paste("Pearson correlation (all data, even (0,0)): ", format(allCor, digits=3, nsmall=3), "", sep='')
              , paste("Pearson correlation omitting all (0,0) points: ", format(oneCanBeZeroCor, digits=3, nsmall=3), "", sep='')
              , paste("Pearson correlation omitting points with ANY zero value: ", format(neitherIsZeroCor, digits=3, nsmall=3), "", sep='')
              , ifelse(log, paste("(Correlation is calculated from the log2(X+1) values, not the original values.)"), '')
              )
            )
}


wget_if_missing <- function(filename, url, verbose=T, comment='') {
  exitcode = 0
  if (file.exists(filename)) {
    if (verbose) { print(paste0("[Not redownloading] existing file: ", filename)); }
    return(TRUE)
  } else {
    exitcode = system(paste0('wget "', url, '" --output-document="', filename, '"'))
  }
  return(file.exists(filename) && (exitcode == 0))
}

system_if_missing <- function(results.vec, ...) {
  # results.vec: result files
  if (all(file.exists(results.vec))) {
    print(paste0("[Not regenerating existing files]: ", paste(results.vec, collapse=", "))); return(TRUE);
  }
  return(system(...))
}

agw_datatable <- function(data, caption="Data Table", options=NULL, ...) {
  require("DT"); # install.packages('DT')
  if (is.null(data) || is.na(data)) { warning("Failure in data table--NA or NULL was passed in!"); stop("Something is messed up, your data to 'agw_datatable' is NULL or NA.") }
  if (is.null(options)) { options=list(pageLength=20, lengthMenu=c(20,50,200,999), autoWidth=TRUE) }
  DT::datatable(data=data, caption=caption, filter=list(position='top',plain=TRUE),  rownames=FALSE, options=options, ...)   #style="bootstrap",
}

agw_read_blast_out6_into_tibble <- function(filename) {      # Reads a blastn '--outfmt 6' output file into a tibble
  require("tibble"); require("readr")  # install.packages("readr") # install.packages("tibble")
  cnames <- c("query_id", "refgenome_id", "perc_ident", "align_length", "n_mismatch", "n_gapopen", "query_start", "query_end", "refgenome_start", "refgenome_end", "evalue", "bit_score")
  btib <- readr::read_tsv(filename, comment="#", na="", col_names=cnames)
  return(btib)
}

agw_read_bed_into_tibble <- function(filename) {
     library("tibble"); library("dplyr") # install.packages("tibble")
     # probably should be read_tsv
     x <- as_tibble(read.table(filename, sep="\t", comment.char='', quote='', header=F, na.strings="", stringsAsFactors=F))
     # Example:  L1PA15	3375	3764	orf	0	-
     #           L1PA15	3820	4179	orf	0	+
     stopifnot(ncol(x) >= 3)
     bedtib <- tibble(chrom    = pull(x[,1])
                   ,  start = as.integer(pull(x[,2]))
                ,    end = as.integer(pull(x[,3]))
               ,  width = abs(as.integer(pull(x[,3])) - as.integer(pull(x[,2]))) # BED is 0-based, so you just subtract (end-start) (see: https://github.com/arq5x/bedtools2/blob/master/docs/content/overview.rst
             ,   name = if(ncol(x)>=4) {            pull(x[,4])  } else { rep("NA", times=nrow(x)) } # name is optional. Default: "NA" (string! Not the R NA special value)
            ,  score = if(ncol(x)>=5) { as.integer(pull(x[,5])) } else { rep(0   , times=nrow(x)) } # score is optional. Default: 0.
            , strand = if(ncol(x)>=6) {            pull(x[,6])  } else { rep('.' , times=nrow(x)) } # strand is optional. '.' means no data
              )
     return(bedtib)
}

agw_convert_blast_tibble_into_bed <- function(x) {
  # Converts one of the "blast" data structures (from 'read_blast_out6_into_tibble') above into a tibble
  return(tibble("chrom"=x$refgenome_id
             , "start"  = x$refgenome_start
             , "end"    = x$refgenome_end
             , "width"  = abs(x$refgenome_end - x$refgenome_start)+1 # looks like 1-100 means "all 100 bp", so we add one here. UNLIKE "read_bed_into_tibble!"
             , "name"   = x$query_id
             , "score"  = rep(1000, times=nrow(x))
             , "strand" = rep('.' , times=nrow(x))
             ))
}

agw_bed_plot <- function(bed_tibble, annot_tib=NULL, genome_len_table=NULL, title="", display_names=TRUE, left_label_space_in_bp=1000) {
  library("ggplot2") # install.packages("ggplot2")
  BAND_HEIGHT_IN_UNITS     = 1.0    # larger regions for each bar
  BAND_FRACTION_SPACE_USED = 0.50   # larger number = fatter bars
  ##MAX_Y_ADJUSTMENT_FOR_ANNOTATION = 0.1 # move the ORFs around a little
  #INTRON_FRACTION_OF_EXON  = 0.85
  LABEL_SPACING_ON_LEFT_IN_BASE_UNITS  = 50
  LABEL_SPACING_ON_RIGHT_IN_BASE_UNITS = LABEL_SPACING_ON_LEFT_IN_BASE_UNITS
  LABEL_MAX_LEFT_EXTENT_IN_BASES       = left_label_space_in_bp # this is very specific to the exact size of the plot unfortunately, so it will need to be auto-set some other way if we aren't looking at repetitive elements

  # The bed_tibble shows the elements you want to plot onto the chromosomes (given in the 'chrom' field of the bed tibble). This is used to visualize BLAST hits.
  # The genome_len_table, if present, shows how long the chromosomes are.
  # The annot_tib, if present, shows the features on your chromosomes.
  bed2 = as_tibble(merge(bed_tibble, genome_len_table, by="chrom", all=F, sort=F))
  bed2[["chrom_idx"]]      = as.integer(as.factor(bed2[["chrom"]]))
  bed2[["y1"]]             = bed2[["chrom_idx"]] * BAND_HEIGHT_IN_UNITS
  bed2[["y2"]]             = bed2[["y1"]] + BAND_FRACTION_SPACE_USED*BAND_HEIGHT_IN_UNITS #*ifelse(is.null(annot_tib), 1.0, INTRON_FRACTION_OF_EXON)
  bed2[["name_truncated"]] = agw_truncate_text(bed2[["name"]])

  if (!is.null(annot_tib)) {
    annot2 <- as_tibble(merge(annot_tib, bed2, by="chrom", all=F, sort=F, suffixes=c("", ".bed")))
    annot2[["name_truncated"]] = agw_truncate_text(annot2[["name"]])
    #annot2[["y2"]] = annot2[["y1"]] + BAND_FRACTION_SPACE_USED*BAND_HEIGHT_IN_UNITS #*(annotation_y_adjustments)
    orf_add_gene_features = geom_rect(data=annot2, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2), color=NA, fill="black", alpha=0.15) # ORF ALPHA
  } else {
    orf_add_gene_features = NULL
  }
  
  bigchroms = data.frame(chrom=bed2$chrom, x1=0, x2=bed2[["chromsize"]], y1=bed2$y1, y2=bed2$y2)[ !duplicated(bed2$y1), ] # the FULL "chromosomes"
  DEFAULT_NAME_COLOR = "#444444" # dark gray
  element_colors = rep(DEFAULT_NAME_COLOR, times=length(bigchroms$chrom))
  element_colors[grepl("^Alu" , bigchroms$chrom, perl=T)  ] = "#00FF00" # green
  # ddd = data.frame(x1=c(10,40,60),y1=c(10,50, 10),x2=c(30,45, 70),y2=c(20,59, 20),wid=c(5,10, 15), widfil=c("ZZred","AAblue","CCorange"), stringsAsFactors=F);   z = ggplot() + geom_rect(data=ddd, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=widfil)) + scale_fill_manual(values=c("red", "blue", "orange")[order(ddd$widfil)] ); print(z)
  text_data = bigchroms
  
  UNSET_COLOR          = "#FF00FF"
  WIDTH_90_PLUS_COLOR  = "#FF0000" # bright red
  WIDTH_80_90_COLOR    = "#EEBB00" # orange-yellow
  WIDTH_UNDER_80_COLOR = "#6778FF" #blue-ish
  #elem_color_levels = c(UNSET_COLOR, WIDTH_90_PLUS_COLOR, WIDTH_80_90_COLOR, WIDTH_UNDER_80_COLOR)
  fillcats=c("UNSET", "90+ bp", "80-90 bp", "<80 bp")
  fillcolors=c(UNSET_COLOR, WIDTH_90_PLUS_COLOR, WIDTH_80_90_COLOR, WIDTH_UNDER_80_COLOR)
  names(fillcolors) <- fillcats
  
  bed2[["fillcategories"]] <- rep(fillcats[1], times=nrow(bed2))
  bed2[["fillcategories"]][ bed2[["width"]] >= 90                        ]   <- fillcats[2]
  bed2[["fillcategories"]][ bed2[["width"]] >= 80 & bed2[["width"]] < 90 ]   <- fillcats[3]
  bed2[["fillcategories"]][                         bed2[["width"]] < 80 ]   <- fillcats[4] # this is SUPER order dependent!!! ggplot is crazy sometimes
  #browser()
  # y1 is the BOTTOM, y2 is the TOP
  ggplot() + 
    scale_x_continuous(name="Position within the element", limits=c(0-LABEL_MAX_LEFT_EXTENT_IN_BASES, NA)) + 
    scale_y_continuous(name="", breaks=NULL) +
    geom_rect(data=bigchroms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0.3) +
    orf_add_gene_features +
    geom_rect(data=bed2, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=fillcategories), alpha=1.0) +
    scale_fill_manual("Probe size (bp)", values=fillcolors) +
    geom_rect(data=bigchroms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0.0) +
    geom_text(data=bed2, aes(x=start, y=y2, label=name_truncated), hjust=0.0, vjust=-0.5, size=2, angle=0) + # <-- PROBE labels. To rotate: vjust=1.0 and angle=90 and hjust=-0.05
    #geom_text(data=bigchroms, aes(x=x1, y=y1+(y2-y1)/2, label=chrom), size=2, hjust=1.0) +
    geom_text(data=text_data, aes(x=x2+LABEL_SPACING_ON_RIGHT_IN_BASE_UNITS, y=y1+(y2-y1)/2, label=paste0(x2, " bp")), color="blue", size=3, hjust=0.0, angle=0) +
    geom_label(data=bigchroms, aes(x=x1-LABEL_SPACING_ON_LEFT_IN_BASE_UNITS, y=y1+(y2-y1)/2, label=chrom), fill=element_colors, color="white", fontface="bold", hjust=1.0, size=4) +
    ggtitle("Probes (red). Dark regions have annotation.")
}

agw_truncate_text <- function(s, n=20, more_char="…") {
  stopifnot(n >= 2)
  truncated = sapply(s, function(sub_s) {
    sarr = unlist(strsplit(as.character(sub_s), split=''))
    if (length(sarr) <= n) {  return(as.character(sub_s)) # just return the original
    } else {                  return(paste0(paste0(head(sarr, n=(n-1)), collapse=''), more_char))  }
  })
  return(unname(truncated))
}

agw_summarize_blast <- function(blast_result_object, reference_fasta, reference_annot_tib, display_names=TRUE, left_label_space_in_bp=1000) {
  library("Biostrings")
  refgen.ss <- Biostrings::readDNAStringSet(filepath=reference_fasta)
  if (nrow(blast_result_object) <= 0) { browser(); }
  bed <- agw_convert_blast_tibble_into_bed(blast_result_object)
  MAX_LEN_FALLBACK <- max(bed$end)
  uniq_chrom_names <- as.character(unique(bed$chrom)) # or levels, since it's a factor
  chrom_lengths.vec <- sapply(uniq_chrom_names, function(x) {
    if (x %in% names(refgen.ss)) { return(width(refgen.ss[x])) }
    else { return(MAX_LEN_FALLBACK); }    })
  #chrom_lengths2.vec <- ifelse(uniq_chrom_names %in% names(refgen.ss), yes=width(refgen.ss[uniq_chrom_names]), no=MAX_LEN_AS_FALLBACK)
  #stopifnot(all(chrom_lengths.vec == chrom_lengths2.vec))
  genome_dat <- data.frame("chrom"=uniq_chrom_names, "chromsize"=chrom_lengths.vec)
  agw_bed_plot(bed_tibble=bed, annot_tib=reference_annot_tib, genome_len_table=genome_dat, display_names=display_names, left_label_space_in_bp=left_label_space_in_bp)
}


### ===============================================================================
### 
### ===============================================================================







