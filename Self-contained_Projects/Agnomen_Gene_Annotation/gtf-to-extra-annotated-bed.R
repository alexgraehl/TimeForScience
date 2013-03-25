# options(error=recover)

library(IRanges) # from Bioconductor. REQUIRED!

trimWhitespaceBothEnds <- function(string) { gsub("^\\s+|\\s+$", '', string) }

print("Ok...")

if (!exists("hugeData")) {
     print("Loading a huge GTF file!")
     hugeData <- read.table("Homo_sapiens.GRCh37.70.gtf.txt", sep="\t"
                            , stringsAsFactors=F
#                            , stringsAsFactors=TRUE
#                            , colClasses=c("factor","factor","factor" # CHROMOSOME, type (e.g. "retained_intron"), exon/CDS/whatever
#                              , "integer", "integer", "character" # START, STOP, SCORE (usually score is just a '.')
#                              , "factor", "factor", "character") # STRAND, something, character
                            )
     print("Loaded the huge GTF file!")
}

print("Setting stuff up...")

#TESTRANGE <- 12250:12255
TESTRANGE <- 1000:10000
tomato <- hugeData[TESTRANGE,]
#tomato <- hugeData

FREETEXT_COL_IDX <- 9

FREETEXT_DELIM <- ";"

print("All set to go...")

## ================= GET THE FREE TEXT COLUMN, WHICH WE WILL BE PROCESSING ================

# gene_id
# transcript_id
# exon_id
# exon_number
# gene_name

xprint <- function(...) {
     print(paste(..., sep=''))
}

xPrintRow <- function(type, range, rowItem) {
     cat(paste(type, rangeStrTab(range), paste(rowItem, collapse="\t"), "\n", sep="\t"))
}


xPrintManualRow <- function(type, smallCoord, bigCoord, rowItem) {
     # prints the two coordinates and the width.
     # (bigCoord - smallCoord + 1) is the width of this feature!
     cat(paste(type, smallCoord, bigCoord, (bigCoord - smallCoord + 1), paste(rowItem, collapse="\t"), "\n", sep="\t"))
}

xPrintManualRowRaw <- function(type, smallCoord, bigCoord, rowItem) {
     # prints the two coordinates and the width.
     # (bigCoord - smallCoord + 1) is the width of this feature!
     #colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")
     cat(paste(type, smallCoord, bigCoord, (bigCoord - smallCoord + 1)
               , paste(c(rowItem[1], rowItem[2], "NA", "NA", "NA", rowItem[6], rowItem[7], rowItem[8], "NA", "NA", rowItem[11]), collapse="\t")
               #, paste(c(rowItem["CHR"], rowItem["ABOUT"], "NA", "NA", "NA", rowItem["STRAND"], rowItem["GENE"], rowItem["TS"], "NA", "NA", rowItem["COMMON"]), collapse="\t")
               , "\n", sep="\t"))
}

xPrintRowRaw <- function(type, range, rowItem) {
     # Print ths row, but without the type, start, stop, and exon/exnum fields.
     # Used for printing introns and UTRs, which don't have valid data for those files in the input row data.
     #                      1   2        3     4       5      6         7     8     9     10      11
     #colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")
     cat(paste(type, rangeStrTab(range)
               , paste(c(rowItem[1], rowItem[2], "NA", "NA", "NA", rowItem[6], rowItem[7], rowItem[8], "NA", "NA", rowItem[11]), collapse="\t")
               #, paste(c(rowItem["CHR"], rowItem["ABOUT"], "NA", "NA", "NA", rowItem["STRAND"], rowItem["GENE"], rowItem["TS"], "NA", "NA", rowItem["COMMON"]), collapse="\t")
               , "\n", sep="\t"))
}


xPrintRowExceptStartStop <- function(type, range, rowItem) {
     # Print ths row with everything EXCEPT the start and stop
     # Used for printing the split-up exons that are partially coding and partially UTR.
     #colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")
     cat(paste(type, rangeStrTab(range)
               , paste(c(rowItem[1], rowItem[2], rowItem[3], "NA", "NA", rowItem[6], rowItem[7], rowItem[8], rowItem[9], rowItem[10], rowItem[11]), collapse="\t")
               #, paste(c(rowItem["CHR"], rowItem["ABOUT"], rowItem["TYPE"], "NA", "NA", rowItem["STRAND"], rowItem["GENE"], rowItem["TS"], rowItem["EXON"], rowItem["EXNUM"], rowItem["COMMON"]), collapse="\t")
               , "\n", sep="\t"))
}

coordsAreEntirelyBetween <- function(start1, end1, start2, end2) {
     # is the ONE item entirely within the TWO item?
     return(start1 >= start2 && end1 <= end2)
}

featureIsEntirelyBetween <- function(irange1, irange2) {
     # range 1 is ENTIRELY within range2 !
     return((irange1@start >= irange2@start) && ((irange1@start+irange1@width) <= (irange2@start+irange2@width)))
     #return(length(IRanges::findOverlaps(irange1, irange2, type="within")) > 0) # Type: WITHIN only
}

#featureIsPartiallyWithin <- function(irange1, irange2) {
#     return(length(IRanges::findOverlaps(irange1, irange2)) > 0)
#}

coordsAreEntirelyOutside <- function(start1, end1, start2, end2) {
     # is the ONE item entirely outside, with no overlap, the TWO item?
     return(end1 < start2 || start1 > end2)
}

featureIsEntirelyOutside <- function(irange1, irange2) {
     ## Either 11111111 22222222222
     ## or:    2222222222 11111111111
     return(((irange1@start+irange1@width-1) < irange2@start) || ((irange2@start+irange2@width-1) < irange1@start))
     #return(0 == length(IRanges::findOverlaps(irange1, irange2)))
}

rangeStr <- function(irange) {
     # returns the string for a range. Example: [1094-2094] (101) . Note that the width is INCLUSIVE! (e.g. 200 to 300 has 101 elements.)
     if (length(irange) == 1) {
          return(paste("[", irange@start, "-", (irange@start+irange@width-1), " width=", irange@width, "]", sep=''))
     } else if (length(irange) == 0) {
          return("[ NULL RANGE ]")
     } else {
          return("[ MULTIPLE IRANGES ARE NOT SUPPORTED BY rangeStr ]");
     }
}


rangeStrTab <- function(irange) {
     # returns the tab-delimited three-element string for a range. Note that the width is INCLUSIVE! (e.g. 2 to 3 has a width of TWO.)
     if (length(irange) == 1) {
          return(paste(irange@start, (irange@start+irange@width-1), irange@width, sep='\t'))
     } else if (length(irange) == 0) {
          return("NA\tNA\tNA")
     } else {
          return("[ MULTIPLE IRANGES ARE NOT SUPPORTED BY rangeStrTab ]\tERROR\tERROR");
     }
}

devour <- function(key, searchInThisVector) {
     # Looks for the KEY (a perl search string!) in the split-up free text vector
     # example input:  devour("name=", someVector) (where someVector could look like this: [1] city=Boston [2] name=Testguy [3] fish=tuna
     # That would return "Testguy" . It returns the NOT_FOUND_PLACEHOLDER_VALUE if there was no match.
     #stopifnot(is.vector(searchInThisVector))
     res <- grep(key, searchInThisVector, perl=T, ignore.case=T)
     if (0 == length(res)) {
          NOT_FOUND_PLACEHOLDER_VALUE <- ""
          return(NOT_FOUND_PLACEHOLDER_VALUE)
     } else {
          stopifnot(1 == length(res)) # length must be EXACTLY one! -- otherwise it means this string was present multiple times, which confuses us!
          return(trimWhitespaceBothEnds(sub(key, '', searchInThisVector[res]))) # delete the key, then trim whitespace again
     }
}

## handleARow <- function(rowItem) {
##      # "type" index is 3
##      # start index is 4, stop index is 5
##      queryRange <- IRanges(start=as.numeric(rowItem[4]), end=as.numeric(rowItem[5]))  # IRanges::intersect(IRanges(10,20), IRanges(40,51))
##      #print( paste("Query: ", codingRange@start, " and ", (codingRange@start + codingRange@width - 1), " is within the coding sequence.", sep=''))
##      if (featureIsEntirelyBetween(queryRange, codingRange)) {
##           xPrintRow(toupper(rowItem[3]), queryRange, rowItem)
##      } else if (featureIsEntirelyOutside(queryRange, codingRange)) {
##           # Any exons BEFORE the start codon are 5' UTR, and any exons AFTER the stop codon are 3' UTR
##           isFivePrimeUtr = (isPositiveStrand && queryRange@start < codingRange@start) || (!isPositiveStrand && queryRange@start > codingRange@start)               # it's either  HERE >>>>>>>>>>>>>>>>>> (on the positive strand)   or                <<<<<<<<<<<<<<<<<<<<<<< HERE (on the negative strand)
##           isThreePrimeUtr = (!isPositiveStrand && queryRange@start < codingRange@start) || (isPositiveStrand && queryRange@start > codingRange@start)               # it's either       HERE <<<<<<<<<<<<<<<<<<<<<<<< (negative strand) or                           >>>>>>>>>>>>>>>>>> HERE   (positive strand)
##           stopifnot(isFivePrimeUtr + isThreePrimeUtr == 1) # EXACTLY one of the two conditions must be true, or there's something wrong in how I programmed this!
##           utrTypeString <- ifelse(isFivePrimeUtr, yes="UTR_FIVE_PRIME", no="UTR_THREE_PRIME")
##           xPrintRowRaw(utrTypeString, queryRange, rowItem)
##      } else {
##           #print(paste("Partially inside at ", rangeStr(queryRange), " ****", sep=''))
##           queryCodeOnly  <- IRanges::intersect(queryRange, codingRange) # only the part of the query that is ALSO in the coding region
##           queryUtrOnly   <-   IRanges::setdiff(queryRange, codingRange) # subtract out the coding part from the query
##           xPrintRowExceptStartStop("EXON_CODING_ONLY",   queryCodeOnly, rowItem) # only the coding part!
##           xPrintRowExceptStartStop("UTR_SECTION_OF_EXON", queryUtrOnly, rowItem) # only the NON-CODING part
##           xPrintRow(paste(toupper(rowItem[3]), "_INCLUDING_UTR", sep=''), queryRange, rowItem) # the exon with both coding AND non-coding parts
##      }
## }

startTime = date() # save the time when we started

BIGGEST_POSSIBLE_GENOMIC_COORDINATE <- (.Machine$integer.max - 100) # This needs to be greater than any valid genomic coordinate, but not so big as to instantly cause integer overflow. What a mess. Why did I use IRanges anyway!!!
NUM_INTERESTING_ITEMS_IN_RESULT_VEC <- 5

xprint("Now setting up the results matrix... this takes a while")
xprint("### Current time is: ", date())

resMat <- matrix("", nrow=nrow(tomato), ncol=NUM_INTERESTING_ITEMS_IN_RESULT_VEC)
column9 <- strsplit(tomato[, FREETEXT_COL_IDX], FREETEXT_DELIM, perl=T)
#column9 <- lapply(colNineWithWhitespace, trimWhitespaceBothEnds) # apply the trimming function to remove unnecessary whitespace
for (i in seq_along(column9)) {
     resMat[i, 1:NUM_INTERESTING_ITEMS_IN_RESULT_VEC] <- c(devour("gene_id "      , column9[[i]])
                                                           , devour("transcript_id ", column9[[i]])
                                                           , devour("exon_id "      , column9[[i]])
                                                           , devour("exon_number "  , column9[[i]])
                                                           , devour("gene_name "    , column9[[i]]))
}

xprint("Done setting up the results matrix...")
xprint("### Current time is: ", date())

stopifnot(1==2)

INTERESTING_ORIGINAL_COLUMNS <- c(1,2,3,4,5,7) # everything but the free text column and useless score columns
ddd <- cbind( tomato[ , INTERESTING_ORIGINAL_COLUMNS], resMat) # smash them together!
colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")

uniqTS.vec <- as.character(unique(ddd$TS)) # unique transcripts only (as characters)

xprint("### About to redirect all future output to the file alex.test.r.txt...")
#sink(file="alex.test.r.txt", append=FALSE)

numRowsProcessed <- 0
for (ut in uniqTS.vec) {
     b.frame <- ddd[ ddd$TS == ut, , drop=F] # we dont' actually even care if the transcripts are presented in order, although it's nice if they are
     tsFullRange <- IRanges(start=min(b.frame$START, b.frame$STOP), end=max(b.frame$START, b.frame$STOP)) # Start it at the FULL transcript range -- min to max
     exons.frame <- b.frame[ b.frame$TYPE %in% c("exon"), c("START","STOP")] # exons only!
     intronRanges <- setdiff(tsFullRange, IRanges(start=exons.frame$START, end=exons.frame$STOP)) # take the FULL range and then subtract out any exons! What's left is the introns.
     
     nStartCodons = sum(b.frame$TYPE == "start_codon")
     nStopCodons  = sum(b.frame$TYPE ==  "stop_codon")
     
     strand = unique(b.frame$STRAND)
     stopifnot(length(strand) == 1)
     isPositiveStrand <- (strand == "+" || strand == "1" || strand == "+1")

     #print(paste("Strand: ", strand, sep=''))
     #xprint("Transcript <", ut, "> had ", nStartCodons, " start codons (", rep("S", times=nStartCodons),  ") and ", nStopCodons, " stop codons (", rep("E", times=nStopCodons), ").")
     #if (nStartCodons != nStopCodons) {
          #print("### WARNING: UNEQUAL NUMBER OF START & STOP CODONS. This appears to be a common situation, so I guess it is by design?")
     #}

     codingLimitLeft  <- -1 # <-- can be any out-of-bounds invalid genomic coordinate
     codingLimitRight <- (codingLimitLeft-1) # one less than codingLimitLeft, so that IRanges is happy. CANNOT be more than one unit less!
     if (nStartCodons > 0) {
          if (isPositiveStrand) {
               codingLimitLeft  <- as.numeric(b.frame[b.frame$TYPE == "start_codon", "START"])
               codingLimitRight <- ifelse(nStopCodons == 0, yes=BIGGEST_POSSIBLE_GENOMIC_COORDINATE, no=as.numeric(b.frame[b.frame$TYPE == "stop_codon", "STOP"])) # BIGGEST_POSSIBLE_GENOMIC_COORDINATE means "no limit: positive infinity"-ish
          } else if (!isPositiveStrand && nStartCodons > 0) {
               codingLimitRight <- as.numeric(b.frame[b.frame$TYPE == "start_codon", "STOP"])
               codingLimitLeft  <- ifelse(nStopCodons == 0, yes=-1, no=as.numeric(b.frame[b.frame$TYPE == "start_codon", "START"])) # -1 means "no limit: negative infinity"
          }
     }

     if (length(codingLimitLeft) > 1 || length(codingLimitRight) > 1) {
          xprint("### WARNING: For transcript ", as.character(ut), ", the number of start codons (", nStartCodons, ") or number of stop codons (", nStopCodons, ") was greater than 1! We are going to just take the EARLIER start codon and LATEST stop codon (for positive strand; for negative, we do it the opposite way).")
          codingLimitLeft <- min(codingLimitLeft)
          codingLimitRight <- max(codingLimitRight) # just pick the "more extreme" of the start/stop codons
     }
     
     codingRange <- IRanges(start=codingLimitLeft, end=codingLimitRight)
     #xprint("Anything within ", rangeStr(codingRange), " is within the coding sequence.")
     
     for (rrr in seq_along(rownames(b.frame))) {
          rowItem <- b.frame[rrr, ]

          qstart <- rowItem[4] # 4 is the start
          qend   <- rowItem[5] # 5 is the end
          
          #queryRange <- IRanges(start=rowItem$START, end=rowItem$STOP)  # IRanges::intersect(IRanges(10,20), IRanges(40,51))
          #print( paste("Query: ", codingRange@start, " and ", (codingRange@start + codingRange@width - 1), " is within the coding sequence.", sep=''))

          if (coordsAreEntirelyBetween(qstart, qend, codingLimitLeft, codingLimitRight)) {
               xPrintRow(toupper(rowItem$TYPE), queryRange, rowItem)
          } else if (coordsAreEntirelyOutside(qstart, qend, codingLimitLeft, codingLimitRight)) {
               # Any exons BEFORE the start codon are 5' UTR, and any exons AFTER the stop codon are 3' UTR
               isFivePrimeUtr  <- ( isPositiveStrand && qstart < codingLimitLeft) || (!isPositiveStrand && qstart > codingLimitLeft)               # it's either  HERE >>>>>>>>>>>>>>>>>> (on the positive strand)   or                <<<<<<<<<<<<<<<<<<<<<<< HERE (on the negative strand)
               isThreePrimeUtr <- (!isPositiveStrand && qstart < codingLimitLeft) || ( isPositiveStrand && qstart > codingLimitLeft)               # it's either       HERE <<<<<<<<<<<<<<<<<<<<<<<< (negative strand) or                           >>>>>>>>>>>>>>>>>> HERE   (positive strand)
               stopifnot(isFivePrimeUtr + isThreePrimeUtr == 1) # EXACTLY one of the two conditions must be true, or there's something wrong in how I programmed this!
               utrTypeString <- ifelse(isFivePrimeUtr, yes="UTR_FIVE_PRIME", no="UTR_THREE_PRIME")
               xPrintManualRowRaw(utrTypeString, qstart, qend, rowItem)
          } else {
               queryRange <- IRanges(start=rowItem$START, end=rowItem$STOP)  # IRanges::intersect(IRanges(10,20), IRanges(40,51))
               #print(paste("Partially inside at ", rangeStr(queryRange), " ****", sep=''))
               queryCodeOnly  <- IRanges::intersect(queryRange, codingRange) # only the part of the query that is ALSO in the coding region
               queryUtrOnly   <-   IRanges::setdiff(queryRange, codingRange) # subtract out the coding part from the query
               xPrintRowExceptStartStop("EXON_CODING_ONLY",   queryCodeOnly, rowItem) # only the coding part!
               xPrintRowExceptStartStop("UTR_SECTION_OF_EXON", queryUtrOnly, rowItem) # only the NON-CODING part
               xPrintRow(paste(toupper(rowItem$TYPE), "_INCLUDING_BOTH_CODING_AND_UTR", sep=''), queryRange, rowItem) # the exon with both coding AND non-coding parts
          }

          numRowsProcessed <- (numRowsProcessed+1)
          if (numRowsProcessed %% 1000 == 0) {
               flush.console()
               #stopifnot(1 == 2)
          }
     }

     #apply(b.frame, 1, handleARow)
     
     for (i in seq_along(intronRanges)) {
          theIntronRange <- intronRanges[i]
          xPrintRowRaw("INTRON", theIntronRange, b.frame[1,])
     }

}

#for (i in seq_along(column9)) {

# Now let's go through them one TRANSCRIPT at a time



xprint("### Start time was: ", startTime)
xprint("### Current time is: ", date())
xprint("### Processed this many features: ", numRowsProcessed)


sink(file=NULL) # stop redirecting printed output to a file!

xprint("### Check the file alex.test.r.txt for output!")
xprint("### Start time was: ", startTime)
xprint("### Current time is: ", date())
xprint("### Processed this many features: ", numRowsProcessed)
