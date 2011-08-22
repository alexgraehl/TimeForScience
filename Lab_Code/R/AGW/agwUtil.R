
## This file should not *source* any others!


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


## =================================================================
## Pass in the NAME of a variable, not a variable!!
## Returns TRUE if there is some data that isn't NA or NULL
## in this structure or variable.
## Returns FALSE if the structure is length zero, if no
## such variable exists, or if the variable is set to NA or NULL.
## So here are some things returning false:
##  list()     NA     NULL     vector(length=0)
## =================================================================
agwHasContent <- function(vName) {
     stopifnot(typeof(vName) == "character")
     if (!exists(vName)) { return(FALSE); }
     v <- eval(parse(text=vName))
     if (length(v) == 0) { return(FALSE); }
     if (is.null(v)) {     return(FALSE); }
     if (typeof(v) != "environment") { ## <-- hashes should not (and cannot) be checked for NA-ness
          if (length(v) == 1 && is.na(v)) { return(FALSE); }
     }
     return(TRUE);
}



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






agwMatrixFromLists <- function(inList) {
     ### Makes a vector with maybe-different-length lists.
     ### Pads the extra values with NA
     longestRowLength <- max(sapply(inList, length))
     tempVec <- vector()
     for (i in 1:length(inList)) {
          thisVecLen <- length(inList[[i]])
          v <- c(inList[[i]], rep(NA, longestRowLength-thisVecLen))
          tempVec <- append(tempVec, v)
     }
     return( matrix(tempVec, nrow=length(inList), ncol=longestRowLength, byrow=TRUE) )
}



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


## ====================================================
## ====================================================
## Prints with logging-to-file behavior by default
log.red.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#FF0000", log=log)
}

log.yellow.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#FFFF00", log=log)
}

log.green.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#00FF00", log=log)
}

log.cyan.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#00FFFF", log=log)
}

log.magenta.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#FF00FF", log=log)
}

log.blue.agw <- function(..., newline=T, log=T) {
     print.color.agw(..., newline=newline, fg="#0000FF", log=log)
}
## ====================================================
## ====================================================

## ====================================================
## ====================================================
## Prints to the console, default is NOT to log to a file
print.red.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#FF0000", log=log)
}

print.yellow.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#FFFF00", log=log)
}

print.green.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#00FF00", log=log)
}

print.cyan.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#00FFFF", log=log)
}

print.magenta.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#FF00FF", log=log)
}

print.blue.agw <- function(..., newline=T, log=F) {
     print.color.agw(..., newline=newline, fg="#0000FF", log=log)
}

## ====================================================
## ====================================================

print.color.agw <- function(..., newline=T, log=F, fg=NULL, bg=NULL) {
     if (require(xterm256)) {
          ## "require" checks the currently loaded packages and doesn't reload code that is already loaded.
          #warning("Note that requiring xterm256 apparently also breaks the ability to inspect functions: you will get an error about \"renderer$translator\".")
          print.agw( style(paste(..., collapse=NULL, sep=''), fg=fg, bg=bg, check.xterm=FALSE)
                    , newline=newline, log=log)
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

## agwPreparePlot <- function(directory=NULL, file=NULL, filetype="png", pointsize=DEFAULT.POINTSIZE, width=NULL, height=NULL, res=SCREEN.RES.DPI, verbose=TRUE, func=NULL, fullPath=NULL, ...) {
##      # Height and width should be given in PIXELS, even for pdfs. If you want to make your 8.5-by-11 pdf the right
##      # size image-wise, you should call it like:  agwPreparePlot(.. , res=72, width=(8.5*72), height=(11*72))
##      agwFinishPlot() ## <-- finish a plot if there was already one
     
##      completeFilename <- NULL
     
##      if (is.null(file) && is.null(fullPath)) {
##           if (verbose) { cat(paste(">>>> agwPreparePlot: Status message: About to plot to the standard screen.\n")) }
          
##      } else {    ## Ok, the user DOES want output to a file...
##           if (!is.null(fullPath)) {
##                completeFilename <- paste(fullPath, '.', filetype, sep='')
##           } else {
##                stopifnot(!is.null(directory))
##                stopifnot(!is.null(file))
##                completeFilename <- paste(directory, '/', file, '.', filetype, sep='')
##           }
          
##           if (verbose) { cat(paste(">>>> agwPreparePlot: Status message: About to plot to the file <", completeFilename, ">...", "\n", sep='')) }
          
##           if (!is.null(directory) && !(file.exists(directory))) {
##                cat(agwGlue("\n\n\n>>>> agwPreparePlot: ERROR: Quitting, because the directory \"", directory, "\" did not exist!\n>>>>        You will need to manually create it!\n"))
##                stop(agwGlue("ERROR: agwPreparePlot: Quitting early because directory \"", directory, "\" does not exist! You will need to manually create it!"))
##           }
##           ## Note: you probably want to call agwFinishPlot after you make your plot!
##           if ("png" == filetype) {
##                if (is.null(width)) { width <- STANDARD.PLOT.WIDTH; }
##                if (is.null(height)) { height <- STANDARD.PLOT.HEIGHT; }
##                png(filename=completeFilename, pointsize=pointsize
##                    , width=width, height=height, res=res)
##           } else if ("pdf" == filetype) {
##                ## width and height are assumed to have been
##                ## given in PIXELS
##                pdf(file=completeFilename, pointsize=pointsize, width=(width/res), height=(height/res))
##           } else {
##                stop("agwPreparePlot: Sorry, only PNG and PDF filetypes are supported for now!!!")
##           }
##      }
     
##      if (!is.null(func)) {  ## If there is a function to call, then call it!
##           func(...) ;
##      }
##      gvGLOBAL.LAST.PLOT <<- file ;
## }

## agwFinishPlot <- function() {
##      if (length(dev.list()) > 0) {
##           dev.off();
##      }
##      gvGLOBAL.LAST.PLOT <<- NULL
## }



## =================================================================

agwGlue <- function(...) {
     return(paste(...,sep='')); ## Just "paste" with no separator
}
## =================================================================


agwPrint <- function(...) {
     cat(agwGlue(..., "\n"))
}

## =================================================================
agwPerlSubi <- function(...) { ## case-insensitive, SINGLE search (not gsub!)
     return(sub(..., ignore.case=TRUE, perl=TRUE))
}
## =================================================================

## =================================================================
agwPerlSub <- function(...) { ## case-sensitive, SINGLE search (not gsub!)
     return(sub(..., ignore.case=FALSE, perl=TRUE))
}
## =================================================================


## =================================================================
agwPerlGSubi <- function(...) { ## case-insensitive, GLOBAL
     return(gsub(..., ignore.case=TRUE, perl=TRUE))
}
## =================================================================

## =================================================================
agwPerlGSub <- function(...) { ## case-sensitive!
     return(gsub(..., ignore.case=FALSE, perl=TRUE))
}
## =================================================================

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

## =================================================================

agwMatrixFromLists <- function(inList) {
     ### Makes a vector with maybe-different-length lists.
     ### Pads the extra values with NA
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
                       , formatted=TRUE
                       , requireZeroExitCode=FALSE ## Aborts on a non-zero exit code
                       , dryrun=FALSE
                       , wait=TRUE
                       , log=FALSE ## Should we log to the output file?
                       , time=FALSE ## Should we print the time?
                       ) {
     
     theCommand <- paste(..., sep=sep) ## Default is to just mash all the input arguments together with no spaces.
     syscallPrefix <- "System Call: ";
     if (!wait)  { syscallPrefix <- "System (Background): " }
     if (dryrun) { syscallPrefix <- "System call dry run (not executed!): " }

     fgColor <- ifelse(dryrun, "#00FFFF" ## cyan = dry run
                       , "#00FF00") ## green = real run
     
     numPrefixSpaces <- 1 + nchar(syscallPrefix) # <-- how many spaces required to indent the following lines of the combined system call properly?
     prefixWhitespace <- paste(rep(' ', numPrefixSpaces), collapse='') ## <-- a whitespace region to indent things properly

     if (time) { print.color.agw(paste("[", Sys.time(), "] ", "\n", sep=''), fg=fgColor, log=log, newline=F) } ## <-- note the lack of a newline!
     
     if (formatted) {
          formattedCmd = gsub("    ", paste("\n", prefixWhitespace, sep=''), theCommand, perl=TRUE) ## <-- four spaces becomes a newline
          formattedCmd = gsub(" [|] ", paste("\n", prefixWhitespace,"| ", sep=''), formattedCmd, perl=TRUE) # <-- a vertical bar surrounded by spaces becomes a newline, then a vertical bar
          print.color.agw(syscallPrefix, formattedCmd, "\n", fg=fgColor, log=log)
     } else {
          print.color.agw(syscallPrefix, theCommand, "\n", sep='', fg=fgColor, log=log)
     }

     exitCode <- 0
     if (!dryrun) {
          exitCode <- system(theCommand, wait=wait)  ## <-- the unix error/success code
     }

     if (time) { print.color.agw(paste("Done at [", Sys.time(), "]", "\n", sep=''), fg=fgColor, log=log, newline=TRUE) }
     
     if (requireZeroExitCode) {
          assert.agw(exitCode == 0, paste("Uh oh! The exit code for the command was NOT zero (zero indicates a success in most utilities)! Instead, the exit code was: ", exitCode, "\n", sep=''))
     }
     
     return(exitCode)
}




