#
#
#
## =================================================================
## ## ## ## ## ## ## ## COMMON R CONSTANTS (NOT VARIABLES!) ## ## ## ## ## ## ## ##
## ## ##
## ##
##

LTY.DASHED.LINE     = 3  ## "LTY: Line TYpe"
LAS.HORIZONTAL.TEXT = 1  ## Text alignment for margin text
LAS.PERPENDICULAR   = 2
LAS.VERTICAL.TEXT   = 3  ## Text alignment for margin text

POS.BOTTOM = 1; POS.LEFT = 2; POS.ABOVE = 3 ; POS.RIGHT = 4 ;

A.X.AXIS = 1
A.Y.AXIS = 2

BY.ROW = 1; BY.COL = 2 ; BY.BOTH = c(1,2)

WORST.P.VALUE = 1

PCH.X = 4
PCH.DIAMOND = 23  ## The hollow diamong plotting character
PCH.SMALL.CIRCLE = 20
PCH.MEDIUM.CIRCLE = 16
PCH.BIG.CIRCLE = 19

PCH.SOLID.BOX = 15
PCH.BOX = PCH.SOLID.BOX

PCH.HOLLOW.UP.TRIANGLE   = 2
PCH.HOLLOW.DOWN.TRIANGLE = 6
PCH.BORDER.UP.TRIANGLE    = 24
PCH.BORDER.DOWN.TRIANGLE  = 25

HIGH.RES.DPI       <- 150
SCREEN.RES.DPI     <- 72

DEFAULT.POINTSIZE  = 12

STANDARD.PLOT.WIDTH  = 750
STANDARD.PLOT.HEIGHT = 750


##
## ##
## ## ##
## ## ## ## ## ## ## ## COMMON R CONSTANTS (NOT VARIABLES!) ## ## ## ## ## ## ## ##
## =================================================================




## ==============================================
## ## ## ## ## ## ## ## COMMON FUNCTIONS ## ## ## ## ## ## ## ##
## ## ##
## ##
##



## =================================================================

agwSrc <- function(filename) {
     source(filename)
}

agwSrcAndRun <- function(filename) {
     # Loads the source file "filename"
     # Then it executes a function with the same name as the filename, except with "_main()" instead of ".R".

     # If that ..._main() function does not exist, then this function stops.

     # If you just want to load a file, use "source(filename)" or "agwSrc"

     agwSrc(filename)
     filenameOnlyFromPath <- gsub(".*/", "", filename) # strip all the slashes, leaving us only with the filename name
     
     mainFunctionName <- gsub(".R", "_main", filenameOnlyFromPath)

     if (!exists(mainFunctionName)) {
          print("Was trying to find the main function named:")
          print(mainFunctionName)
          stopifnot(exists(mainFunctionName))
     }
     eval(parse(text=paste(mainFunctionName, "()", sep=''))) # <-- actually call the "filename"_main() function!
}


## =================================================================
agwMax <- function(...) {
     return(max(..., na.rm=TRUE))
}

agwMin <- function(...) {
     return(min(..., na.rm=TRUE))
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

## Like the built-in R "grep", but it returns the actual
## elements found, instead of the indices  where the elements
## were found. (So it is like UNIX grep in this regard.)
agwGrep <- function(what, where, ...) {
     z <- grep(what, where, ...)
     return(where[z])
}

## Same as "agwGrep" above, but ignoring case
agwGrepi <- function(what, where, ...) {
     return(agwGrep(what, where, ..., ignore.case=TRUE))
}

## =================================================================


agwPreparePlot <- function(directory=NULL, file=NULL, filetype="png", pointsize=DEFAULT.POINTSIZE, width=NULL, height=NULL, res=SCREEN.RES.DPI, verbose=TRUE, func=NULL, fullPath=NULL, ...) {
     # Height and width should be given in PIXELS, even for pdfs. If you want to make your 8.5-by-11 pdf the right
     # size image-wise, you should call it like:  agwPreparePlot(.. , res=72, width=(8.5*72), height=(11*72))
     agwFinishPlot() ## <-- finish a plot if there was already one
     
     completeFilename <- NULL
     
     if (is.null(file) && is.null(fullPath)) {
          if (verbose) { cat(paste(">>>> agwPreparePlot: Status message: About to plot to the standard screen.\n")) }
          
     } else {    ## Ok, the user DOES want output to a file...
          if (!is.null(fullPath)) {
               completeFilename <- paste(fullPath, '.', filetype, sep='')
          } else {
               stopifnot(!is.null(directory))
               stopifnot(!is.null(file))
               completeFilename <- paste(directory, '/', file, '.', filetype, sep='')
          }
          
          if (verbose) { cat(paste(">>>> agwPreparePlot: Status message: About to plot to the file <", completeFilename, ">...", "\n", sep='')) }
          
          if (!is.null(directory) && !(file.exists(directory))) {
               cat(agwGlue("\n\n\n>>>> agwPreparePlot: ERROR: Quitting, because the directory \"", directory, "\" did not exist!\n>>>>        You will need to manually create it!\n"))
               stop(agwGlue("ERROR: agwPreparePlot: Quitting early because directory \"", directory, "\" does not exist! You will need to manually create it!"))
          }
          ## Note: you probably want to call agwFinishPlot after you make your plot!
          if ("png" == filetype) {
               if (is.null(width)) { width <- STANDARD.PLOT.WIDTH; }
               if (is.null(height)) { height <- STANDARD.PLOT.HEIGHT; }
               png(filename=completeFilename, pointsize=pointsize
                   , width=width, height=height, res=res)
          } else if ("pdf" == filetype) {
               ## width and height are assumed to have been
               ## given in PIXELS
               pdf(file=completeFilename, pointsize=pointsize, width=(width/res), height=(height/res))
          } else {
               stop("agwPreparePlot: Sorry, only PNG and PDF filetypes are supported for now!!!")
          }
     }
     
     if (!is.null(func)) {  ## If there is a function to call, then call it!
          func(...) ;
     }
     gvGLOBAL.LAST.PLOT <<- file ;
}

agwFinishPlot <- function() {
     if (length(dev.list()) > 0) {
          dev.off();
     }
     gvGLOBAL.LAST.PLOT <<- NULL
}

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
     if (max < min) { stop("ERROR in arguments to agwClamp: max is greater than min!"); }
     if (!is.null(min)) { x[x < min] = min }
     if (!is.null(max)) { x[x > max] = max }
     return(x);
}
## =================================================================

## =================================================================
agwReadFileIntoDataFrame <- function(filename, allowRagged=FALSE, row.names=1, header=TRUE) {
     ## Input: a tab-delimited file with a single row header AND a single column header.
     ## Set <row.names> to NULL to force auto-numbering of rows (useful if you don't care about the row names, and want to avoid the "duplicate 'row.names' are not allowed" problem)
     return(read.table(filename, stringsAsFactors=FALSE
                       , header=header
                       , comment.char=""
                       , sep="\t"
                       , row.names=row.names
                       , quote=""
                       , na.strings=c("NA","NaN","ND")
                       , check.names=FALSE
                       , fill=allowRagged));
}
## =================================================================

agwReadFileIntoListOfLists <- function(filename, sep="\t", row.name.column=NA) {
     ## Reads from a file into a list of lists
     ## "row.name.column" lets you pick which column becomes the row name.
     ## Note that I think R requires that row names are unique.
     ## If "row.name.column" is set to NA, then no row names are set.
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

gvDEPENDENCY.NETWORK <- agwNewHash()
agwAddRule <- function(target=NULL, dependsOn=NULL) {
     ## Prereqs can be a vector or a single item.
     ## It should be textual! Everything here has to be also defined
     ## in agwGlobalLoad or it won't work!
     existingDependencies = agwHashGet(gvDEPENDENCY.NETWORK, target)
     if (is.null(existingDependencies)) {
          print("Adding a new rule...")
     } else {
          print("Modifying an old rule...")
     }
     agwHashPut(gvDEPENDENCY.NETWORK, target, union(existingDependencies, dependsOn))
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
     if (is.null(v)) {     return(FALSE); }
     if (typeof(v) != "environment") { ## <-- hashes should not (and cannot) be checked for NA-ness
          if (length(v) == 1 && is.na(v)) { return(FALSE); }
     }
     return(TRUE);
}

## =================================================================

agwFillRect <- function(size=NA, col="#FF000055") {
     # Basically just fill the background of a plot.
     # You should pick "size" such that it's a big number
     # that covers the biggest extent of the graph.
     rect(xleft=-size, xright=size, ybottom=-size, ytop=size, col=col, border=NA)
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

##
## ##
## ## ##
## ## ## ## ## ## ## ## COMMON FUNCTIONS ## ## ## ## ## ## ## ##
## ==============================================


## Name summing: sum(sapply(strsplit("abgcdefg",'')[[1]],function(x){as.integer(charToRaw(x))}))






##sum(sapply(strsplit(colnames(tplot.mat.sorted.by.max)[ncol(tplot.mat.sorted.by.max)],'')[[1]],function(x){as.integer(charToRaw(x))}))
##strAmt <- sum(sapply(strsplit("abgcdefgg",'')[[1]],function(x){as.integer(charToRaw(x))}))
##sapply(strAmt,function(a){ rgb((a*37) %% 255, (a*13) %% 255, (a*17) %% 255, maxColorValue=255) })
##rgb((strAmt*37) %% 255, (strAmt*13) %% 255, (strAmt*17) %% 255, maxColorValue=255)



