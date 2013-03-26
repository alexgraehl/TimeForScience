# options(error=recover)

library(IRanges) # from Bioconductor. REQUIRED!

OUTFILE <- "alex.test.r.txt"




xprint <- function(...) { cat(paste(..., sep=''), "\n") }
xPrintGeneric    <- function(type, start, end, rowStuff) { cat(paste(type, start, end, (end-start+1), paste(rowStuff, collapse="\t"), sep="\t"), "\n") }
xPrintExonSubset <- function(type, featStart, featEnd, rowItem) {
     #                      1   2        3     4       5      6         7     8     9     10      11
     #colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")
     xPrintGeneric(type=type, start=featStart, end=featEnd, rowStuff=c(rowItem[1:3], "NA", "NA", rowItem[6:11]))
}
xPrintFullExon <- function(type, featStart, featEnd, rowItem) { xPrintGeneric(type, featStart, featEnd, rowItem) }
xPrintUTR      <- function(type, featStart, featEnd, rowItem) { xPrintGeneric(type, featStart, featEnd, c(rowItem[1:2], "NA", "NA", "NA", rowItem[6:8], "NA", "NA", rowItem[11])) }
xPrintIntron   <- function(...) { xPrintUTR(...) } # same for now

coordsAreEntirelyBetween <- function(start1, end1, start2, end2) { return(start1 >= start2 && end1 <= end2) } # is the "start1 & end1" entirely between the "start2 & end2". Note that start is ALWAYS smaller than end in these conventions---no strandedness strangeness is going on here.
coordsAreEntirelyOutside <- function(start1, end1, start2, end2) { return(end1 < start2 || start1 > end2) } # is the "start1 & end1" entirely OUTSIDE the "start2 & end2". Note that start is ALWAYS smaller than end in these conventions---no strandedness strangeness is going on here.

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

devour <- function(key, searchInThisVector) {
     # Looks for the KEY (a perl-style regexp search expression!) in the split-up free text vector
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


BIGGEST_POSSIBLE_GENOMIC_COORDINATE <- (.Machine$integer.max - 100) # This needs to be greater than any valid genomic coordinate, but not so big as to instantly cause integer overflow. What a mess. Why did I use IRanges anyway!!!
NUM_INTERESTING_ITEMS_IN_RESULT_VEC <- 5

if (!exists("hugeData")) {
     print("Loading a huge GTF file!")
     hugeData <- read.table("Homo_sapiens.GRCh37.70.gtf.txt", sep="\t"
                            , stringsAsFactors=F  # Factors just seems to make things SLOWER amazingly
                            , colClasses=c("character","character","character" # CHROMOSOME, type (e.g. "retained_intron"), exon/CDS/whatever
                              , "integer", "integer", "character" # START, STOP, SCORE (usually score is just a '.')
                              , "character", "character", "character") # STRAND, something, character
                            , nrows=50000)
     warning("DEBUG: ONLY LOADED 50,000 rows from hugeData!!!!!!!!!!!!!!!")
     print("Loaded the huge GTF file!")
}

warning("DEBUG: ONLY LOADED 50,000 rows from hugeData!!!!!!!!!!!!!!!")

print("Setting stuff up...")

xprint("By the way, there are a total of ", nrow(hugeData), " rows that we will probably have to process in all.")
#TESTRANGE <- 12250:12255
TESTRANGE <- 1000:1200
tomato <- hugeData[TESTRANGE,]
#tomato <- hugeData

FREETEXT_COL_IDX <- 9

FREETEXT_DELIM <- ";"

print("All set to go...")

startTime = date() # save the time when we started

xprint("Now setting up the results matrix... this takes a while")
xprint("### Current time is: ", date())

resMat <- matrix("", nrow=nrow(tomato), ncol=NUM_INTERESTING_ITEMS_IN_RESULT_VEC)
column9 <- strsplit(tomato[, FREETEXT_COL_IDX], FREETEXT_DELIM, perl=T)
for (i in seq_along(column9)) {
     resMat[i, 1:NUM_INTERESTING_ITEMS_IN_RESULT_VEC] <- c(devour("gene_id "      , column9[[i]])
                                                           , devour("transcript_id ", column9[[i]])
                                                           , devour("exon_id "      , column9[[i]])
                                                           , devour("exon_number "  , column9[[i]])
                                                           , devour("gene_name "    , column9[[i]]))
}

xprint("Done setting up the results matrix...")
xprint("### Current time is: ", date())
xprint("### Mashing the two matrices together now this takes a while... ", date())
INTERESTING_ORIGINAL_COLUMNS <- c(1,2,3,4,5,7) # <-- this means "everything but the free text column and useless score columns!"
ddd           <- cbind( tomato[ , INTERESTING_ORIGINAL_COLUMNS], resMat, stringsAsFactors=FALSE) # smash them together!
colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")

xprint("### Done mashing the two matrices together... ", date())
xprint("### About to redirect all future output to the file ", OUTFILE)
#sink(file=OUTFILE, append=FALSE)
xprint("### Started this file generation at ", date())

numRowsProcessed <- 0
uniqTS.vec <- as.character(unique(ddd$TS)) # unique transcripts only (as characters)
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
          xprint("### WARNING: For transcript ", as.character(ut), ", the number of start codons (", nStartCodons, ") or number of stop codons (", nStopCodons, ") was greater than 1! We are going to just going to use the start/stop pair that covers the most area.")
          codingLimitLeft <- min(codingLimitLeft)
          codingLimitRight <- max(codingLimitRight) # just pick the "more extreme" of the start/stop codons
     }
     
     codingRange <- IRanges(start=codingLimitLeft, end=codingLimitRight)
     for (rrr in seq_along(rownames(b.frame))) {
          #rowItem <- b.frame[rrr,] # list
          if ("CDS" == b.frame[rrr,3] ) { # column 3 is the TYPE
               next; # skip the "CDS" entries---they appear to be totally redundant with the exons, and contribute no useful information... maybe
          }
          qstart <- b.frame[rrr,4] # 4 is the start index. Is it faster to address it as a '4' instead of "START" ? I don't actually know. qstart = QUERY START.
          qend   <- b.frame[rrr,5] # 5 is the end index. qend = QUERY END
          if (coordsAreEntirelyBetween(qstart, qend, codingLimitLeft, codingLimitRight)) {
               xPrintFullExon(toupper(b.frame[rrr,"TYPE"]), qstart, qend, b.frame[rrr,])
          } else if (coordsAreEntirelyOutside(qstart, qend, codingLimitLeft, codingLimitRight)) {
               # Any exons BEFORE the start codon are 5' UTR, and any exons AFTER the stop codon are 3' UTR
               isFivePrimeUtr  <- ( isPositiveStrand && qstart < codingLimitLeft) || (!isPositiveStrand && qstart > codingLimitLeft)               # it's either  HERE >>>>>>>>>>>>>>>>>> (on the positive strand)   or                <<<<<<<<<<<<<<<<<<<<<<< HERE (on the negative strand)
               isThreePrimeUtr <- (!isPositiveStrand && qstart < codingLimitLeft) || ( isPositiveStrand && qstart > codingLimitLeft)               # it's either       HERE <<<<<<<<<<<<<<<<<<<<<<<< (negative strand) or                           >>>>>>>>>>>>>>>>>> HERE   (positive strand)
               stopifnot(isFivePrimeUtr + isThreePrimeUtr == 1) # EXACTLY one of the two conditions must be true, or there's something wrong in how I programmed this!
               utrTypeString <- ifelse(isFivePrimeUtr, yes="UTR_5_PRIME", no="UTR_3_PRIME")
               xPrintUTR(utrTypeString, qstart, qend, b.frame[rrr,])
          } else {
               queryRange <- IRanges(start=b.frame[rrr,"START"], end=b.frame[rrr,"STOP"])  # IRanges::intersect(IRanges(10,20), IRanges(40,51))
               #print(paste("Partially inside at ", rangeStr(queryRange), " ****", sep=''))
               queryCodeOnly  <- IRanges::intersect(queryRange, codingRange) # only the part of the query that is ALSO in the coding region
               queryUtrOnly   <-   IRanges::setdiff(queryRange, codingRange) # subtract out the coding part from the query
               xPrintExonSubset("EXON--CODING_SECTION_ONLY", queryCodeOnly@start, IRanges::end(queryCodeOnly), b.frame[rrr,]) # only the coding part!
               xPrintExonSubset("EXON--UTR_SECTION_ONLY", queryUtrOnly@start, IRanges::end(queryUtrOnly), b.frame[rrr,]) # only the NON-CODING part
               xPrintFullExon(paste(toupper(b.frame[rrr,"TYPE"]), "--COMPLETE--BOTH_CODING_AND_UTR", sep=''), queryRange@start, IRanges::end(queryRange), b.frame[rrr,]) # the exon with both coding AND non-coding parts
          }

          numRowsProcessed <- (numRowsProcessed+1)
          if (numRowsProcessed %% 1000 == 0) {
               xprint("#### DEBUG: processed a total of ", numRowsProcessed, " rows by time ", date())
               flush.console()
          }
     }
     for (i in seq_along(intronRanges)) {
          theIntronRange <- intronRanges[i]
          xPrintIntron("INTRON", theIntronRange@start, end(theIntronRange), b.frame[1,])
     }
}

xprint("### Start time was: ", startTime)
xprint("### Current time is: ", date())
xprint("### Processed this many features: ", numRowsProcessed)

sink(file=NULL) # stop redirecting printed output to a file!

## ========== PRINT TO CONSOLE (after we remove the 'sink') ===============
xprint("### Check the file ", OUTFILE, " for output!")
xprint("### Start time was: ", startTime)
xprint("### Current time is: ", date())
xprint("### Processed this many features: ", numRowsProcessed)




















## rangeStrTab <- function(irange) {
##      # returns the tab-delimited three-element string for a range. Note that the width is INCLUSIVE! (e.g. 2 to 3 has a width of TWO.)
##      if (length(irange) == 1) {
##           return(paste(irange@start, (irange@start+irange@width-1), irange@width, sep='\t'))
##      } else if (length(irange) == 0) {
##           return("NA\tNA\tNA")
##      } else {
##           return("[ MULTIPLE IRANGES ARE NOT SUPPORTED BY rangeStrTab ]\tERROR\tERROR");
##      }
## }

##trimWhitespaceBothEnds <- function(string) { gsub("^\\s+|\\s+$", '', string) }



#featureIsEntirelyBetween <- function(irange1, irange2) {
     # range 1 is ENTIRELY within range2 !
#     return((irange1@start >= irange2@start) && ((irange1@start+irange1@width) <= (irange2@start+irange2@width)))
     #return(length(IRanges::findOverlaps(irange1, irange2, type="within")) > 0) # Type: WITHIN only
#}

#featureIsPartiallyWithin <- function(irange1, irange2) {
#     return(length(IRanges::findOverlaps(irange1, irange2)) > 0)
#}

#featureIsEntirelyOutside <- function(irange1, irange2) {
#     ## Either 11111111 22222222222
#     ## or:    2222222222 11111111111
#     return(((irange1@start+irange1@width-1) < irange2@start) || ((irange2@start+irange2@width-1) < irange1@start))
#     #return(0 == length(IRanges::findOverlaps(irange1, irange2)))
#}
