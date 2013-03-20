# options(error=recover)

library(IRanges) # from Bioconductor. REQUIRED!

trimWhitespaceBothEnds <- function(string) { gsub("^\\s+|\\s+$", '', string) }

print("Ok...")

if (!exists("hugeData")) {
     print("Loading a huge GTF file!")
     hugeData <- read.table("Homo_sapiens.GRCh37.70.gtf.txt", sep="\t", stringsAsFactors=FALSE)
     print("Loaded the huge GTF file!")
}

#TESTRANGE <- 12250:12255
TESTRANGE <- 12250:12355
tomato <- hugeData[TESTRANGE,]


FREETEXT_COL_IDX <- 9

FREETEXT_DELIM <- ";"

## ================= GET THE FREE TEXT COLUMN, WHICH WE WILL BE PROCESSING ================
column9 <- strsplit(tomato[, FREETEXT_COL_IDX], FREETEXT_DELIM, perl=T)
#column9 <- lapply(colNineWithWhitespace, trimWhitespaceBothEnds) # apply the trimming function to remove unnecessary whitespace

# gene_id
# transcript_id
# exon_id
# exon_number
# gene_name

xprint <- function(...) {
     print(paste(..., sep=''))
}

xPrintRow <- function(type, range, rowItem) {
     cat(paste(type, rangeStr(range), paste(rowItem, collapse="\t"), "\n", sep="\t"))
}

xPrintRowRaw <- function(type, range, rowItem) {
     # Print ths row, but without the type, start, stop, and exon/exnum fields.
     # Used for printing introns and UTRs, which don't have valid data for those files in the input row data.
     
     #colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")
     cat(paste(type, rangeStr(range)
               , paste(c(rowItem$CHR, rowItem$ABOUT, "NA", "NA", "NA", rowItem$STRAND, rowItem$GENE, rowItem$TS, "NA", "NA", rowItem$COMMON), collapse="\t")
               , "\n", sep="\t"))
}

featureIsEntirelyBetween <- function(irange1, irange2) {
     return(length(IRanges::findOverlaps(irange1, irange2, type="within")) > 0) # Type: WITHIN only
}

featureIsPartiallyWithin <- function(irange1, irange2) {
     return(length(IRanges::findOverlaps(irange1, irange2)) > 0)
}

featureIsEntirelyOutside <- function(irange1, irange2) {
     return(0 == length(IRanges::findOverlaps(irange1, irange2)))
}

rangeStr <- function(irange) {
     # returns the string for a range. Example: [1094-2094] (101) . Note that the width is INCLUSIVE! (e.g. 200 to 300 has 101 elements.)
     if (length(irange) == 0) {
          return("[ NULL RANGE ]")
     } else if (length(irange) == 1) {
          return(paste("[", irange@start, "-", (irange@start+irange@width-1), " width=", irange@width, "]", sep=''))
     } else {
          return("[ MULTIPLE IRANGES ARE NOT SUPPORTED BY rangeStr ]");
     }
}

devour <- function(key, searchInThisVector) {
     # Looks for the KEY (a perl search string!) in the split-up free text vector
     # example input:  devour("name=", someVector) (where someVector could look like this: [1] city=Boston [2] name=Testguy [3] fish=tuna
     # That would return "Testguy" . It returns the NOT_FOUND_PLACEHOLDER_VALUE if there was no match.
     stopifnot(is.vector(searchInThisVector))
     NOT_FOUND_PLACEHOLDER_VALUE <- ""
     res <- grep(key, searchInThisVector, perl=T, ignore.case=T)
     if (0 == length(res)) {
          return(NOT_FOUND_PLACEHOLDER_VALUE)
     } else {
          stopifnot(1 == length(res)) # length must be EXACTLY one! -- otherwise it means this string was present multiple times, which confuses us!
          return(trimWhitespaceBothEnds(sub(key, '', searchInThisVector[res]))) # delete the key, then trim whitespace again
     }
}

NUM_INTERESTING_ITEMS_IN_RESULT_VEC <- 5
resMat <- matrix("", nrow=nrow(tomato), ncol=NUM_INTERESTING_ITEMS_IN_RESULT_VEC)
for (i in seq_along(column9)) {
     resMat[i, 1:NUM_INTERESTING_ITEMS_IN_RESULT_VEC] <- c("gene_id"        = devour("gene_id "      , column9[[i]])
                                                           , "transcript_id"= devour("transcript_id ", column9[[i]])
                                                           , "exon_id"      = devour("exon_id "      , column9[[i]])
                                                           , "exon_number"  = devour("exon_number "  , column9[[i]])
                                                           , "gene_name"    = devour("gene_name "    , column9[[i]]))
     #print(resMat[i,])
}

INTERESTING_ORIGINAL_COLUMNS <- c(1,2,3,4,5,7) # everything but the free text column and useless score columns
ddd <- cbind( tomato[ , INTERESTING_ORIGINAL_COLUMNS], resMat) # smash them together!
colnames(ddd) <- c("CHR","ABOUT","TYPE","START","STOP","STRAND","GENE","TS","EXON","EXNUM","COMMON")

uniqTS.vec <- as.character(unique(ddd$TS)) # unique transcripts only (as characters)
firstUniqueTranscript <- uniqTS.vec[1]

for (ut in uniqTS.vec) {
     a.frame <- ddd[ ddd$TS == ut, , drop=F]
     b.frame <- a.frame[ order(a.frame$START), ] # re-order the data frame by start positions of each exon/element! It should already be in that order, but just in case...
     ec.frame <- b.frame[ b.frame$TYPE %in% c("exon", "CDS"), c("START","STOP")] # exons and CDs
     #print("Got the 'ec.frame', which contains the start and stop locations for all EXON and CDS elements in this transcript.")

     tsFullRange <- IRanges(start=min(b.frame$START, b.frame$STOP), end=max(b.frame$START, b.frame$STOP)) # Start it at the FULL transcript range -- min to max
     intronRanges <- tsFullRange
     for (i in seq_along(rownames(ec.frame))) {
          # All regions that are 1) NOT in exons
          #              and are 2) between the first and last exons
          #                         are introns.
          intronRanges <- setdiff(intronRanges, IRanges(start=ec.frame[i, "START"], end=ec.frame[i, "STOP"])) # Now SUBTRACT OUT each "exon or CDS" IRange, leaving only things that were NOT in an exon and were NOT in a CDS
     }
     
     nStartCodons = sum(b.frame$TYPE == "start_codon")
     nStopCodons  = sum(b.frame$TYPE ==  "stop_codon")
     
     strand = unique(b.frame$STRAND)
     stopifnot(length(strand) == 1)
     isPositiveStrand <- (strand == "+" || strand == "1" || strand == "+1")

     #print(paste("Strand: ", strand, sep=''))
     xprint("Transcript <", ut, "> had ", nStartCodons, " start codons (", rep("S", times=nStartCodons),  ") and ", nStopCodons, " stop codons (", rep("E", times=nStopCodons), ").")

     if (nStartCodons != nStopCodons) {
          print("WARNING: UNEQUAL NUMBER OF START & STOP CODONS. This appears to be a common situation, so I guess it is by design?")
          #stopifnot(nStartCodons == nStopCodons)
     }
     
     
     BIGGEST_POSSIBLE_GENOMIC_COORDINATE <- (.Machine$integer.max - 100) # This needs to be greater than any valid genomic coordinate, but not so big as to instantly cause integer overflow. What a mess. Why did I use IRanges anyway!!!
     
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

     codingRange <- IRanges(start=codingLimitLeft, end=codingLimitRight)
     xprint("Anything within ", rangeStr(codingRange), " is within the coding sequence.")
     
     for (rrr in seq_along(rownames(b.frame))) {
          rowItem <- b.frame[rrr, ]
          queryRange <- IRanges(start=rowItem$START, end=rowItem$STOP)  # IRanges::intersect(IRanges(10,20), IRanges(40,51))
          #print( paste("Query: ", codingRange@start, " and ", (codingRange@start + codingRange@width - 1), " is within the coding sequence.", sep=''))
          if (featureIsEntirelyBetween(queryRange, codingRange)) {
               xPrintRow(rowItem$TYPEE, queryRange, rowItem)
          } else if (featureIsEntirelyOutside(queryRange, codingRange)) {
               # Any exons BEFORE the start codon are 5' UTR, and any exons AFTER the stop codon are 3' UTR
               # Note that we just check the relative positions of the queryRange start and codingRange start; we don't actually bother to figure out the widths, because we don't care---
               # we ALREADY know that the queryRange is totally outside of the codingRange
               isFivePrimeUtr = (isPositiveStrand && queryRange@start < codingRange@start) || (!isPositiveStrand && queryRange@start > codingRange@start)
               # it's either  HERE >>>>>>>>>>>>>>>>>> (on the positive strand)
               # or                <<<<<<<<<<<<<<<<<<<<<<< HERE (on the negative strand)
               isThreePrimeUtr = (!isPositiveStrand && queryRange@start < codingRange@start) || (isPositiveStrand && queryRange@start > codingRange@start)
               # it's either       HERE <<<<<<<<<<<<<<<<<<<<<<<< (negative strand)
               # or                           >>>>>>>>>>>>>>>>>> HERE   (positive strand)
               stopifnot(isFivePrimeUtr + isThreePrimeUtr == 1) # EXACTLY one of the two conditions must be true, or there's something wrong in how I programmed this!
               utrTypeString <- ifelse(isFivePrimeUtr, yes="UTR_FIVE_PRIME", no="UTR_THREE_PRIME")
               xPrintRowRaw(utrTypeString, queryRange, rowItem)
          } else {
               print(paste("Partially inside at ", rangeStr(queryRange), " ****", sep=''))
               theClassification <- "MIXED CODING/UTR"
               queryCodeOnly  <- IRanges::intersect(queryRange, codingRange) # only the part of the query that is ALSO in the coding region
               queryUtrOnly   <-   IRanges::setdiff(queryRange, codingRange) # subtract out the coding part from the query
               xPrintRowRaw("Coding part of that:", queryCodeOnly, rowItem)
               xPrintRowRaw("UTR part of that: ", queryUtrOnly, rowItem)
               #newCodeRow <- rowItem ; newCodeRow$START <- queryCodeOnly@start; newCodeRow$STOP <- (queryCodeOnly@start+queryCodeOnly@width-1)
               #newUtrRow  <- rowItem ; newUtrRow$START  <- queryUtrOnly@start;  newUtrRow$STOP  <- (queryUtrOnly@start+queryUtrOnly@width-1)
          }
     }

     for (i in seq_along(intronRanges)) {
          theIntronRange <- intronRanges[i]
          xPrintRowRaw("INTRON", theIntronRange, b.frame[1,])
     }
     
}

#for (i in seq_along(column9)) {

# Now let's go through them one TRANSCRIPT at a time











