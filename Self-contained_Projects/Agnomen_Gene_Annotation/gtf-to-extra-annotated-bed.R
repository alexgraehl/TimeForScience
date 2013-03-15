# options(error=recover)

library(IRanges) # from Bioconductor. REQUIRED!

trimWhitespaceBothEnds <- function(string) { gsub("^\\s+|\\s+$", '', string) }

print("Ok...")

if (!exists("hugeData")) {
     print("Loading a huge GTF file!")
     hugeData <- read.table("Homo_sapiens.GRCh37.70.gtf.txt", sep="\t", stringsAsFactors=FALSE)
     print("Loaded the huge GTF file!")
}

TESTRANGE <- 12000:13000
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
     return(paste("[", irange@start, "-", (irange@start+irange@width-1), " width=", irange@width, "]", sep=''))
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
     b.frame <- a.frame[ order(a.frame$START), ] # re-order the data frame by start positions of each exon/element!
     nStartCodons = sum(b.frame$TYPE == "start_codon")
     nStopCodons  = sum(b.frame$TYPE == "stop_codon")

     strand = unique(b.frame$STRAND)
     isPositiveStrand = (strand == "+" || strand == "1" || strand == "+1")
     stopifnot(length(strand) == 1)

     print(paste("Strand: ", strand, sep=''))
     
     print( paste("Transcript <", ut, "> had ", nStartCodons, " start codons (", rep("S", times=nStartCodons),  ") and ", nStopCodons, " stop codons (", rep("E", times=nStopCodons), ").", sep='') )

     if (nStartCodons != nStopCodons) {
          print("WARNING: UNEQUAL NUMBER OF START & STOP CODONS. This appears to be a common situation, so I guess it is by design?")
          #stopifnot(nStartCodons == nStopCodons)
     }
     
     # Any exons BEFORE the start codon are 5' UTR
     
     # Any exons AFTER the stop codon are 3' UTR

     # All regions that are 1) NOT in exons
     #              and are 2) between the first and last exons
     #                         are introns.
     
     #lowerLimit <- min( b.frame$START, b.frame$END ) # lowest-numbered base pair of ANY type. Disregards strand.
     #upperLimit <- max( b.frame$START, b.frame$END ) # highest-numbered base pair of ANY type. Disregards strand.

     BIGGEST_POSSIBLE_GENOMIC_COORDINATE <- (.Machine$integer.max - 100) # This needs to be greater than any valid genomic coordinate, but not so big as to instantly cause integer overflow. What a mess. Why did I use IRanges anyway!!!
     
     codingLimitLeft  <- -1 # any invalid genomic coordinate
     codingLimitRight <- (codingLimitLeft-1) # one less than codingLimitLeft, so that IRanges is happy. CANNOT be more than one unit less!
     if (isPositiveStrand && nStartCodons > 0) {
          codingLimitLeft  = as.numeric(b.frame[b.frame$TYPE == "start_codon", "START"])
          codingLimitRight = ifelse(nStopCodons == 0, yes=BIGGEST_POSSIBLE_GENOMIC_COORDINATE, no=as.numeric(b.frame[b.frame$TYPE == "stop_codon", "STOP"]))
     } else if (!isPositiveStrand && nStartCodons > 0) {
          codingLimitRight = as.numeric(b.frame[b.frame$TYPE == "start_codon", "STOP"])
          codingLimitLeft  = ifelse(nStopCodons == 0, yes=-1, no=as.numeric(b.frame[b.frame$TYPE == "start_codon", "START"]))
     }

     codingRange <- IRanges(start=codingLimitLeft, end=codingLimitRight)
     print( paste("Anything within ", rangeStr(codingRange), " is within the coding sequence.", sep=''))
     
     for (rrr in seq_along(rownames(b.frame))) {
          rowItem <- b.frame[rrr, ]
          queryRange <- IRanges(start=rowItem$START, end=rowItem$STOP)  # IRanges::intersect(IRanges(10,20), IRanges(40,51))

          #print( paste("Query: ", codingRange@start, " and ", (codingRange@start + codingRange@width - 1), " is within the coding sequence.", sep=''))
          if (featureIsEntirelyBetween(queryRange, codingRange)) {
               print(paste("Feature at ", rangeStr(queryRange), " is entirely within the coding limits!", sep=''))
          } else if (featureIsEntirelyOutside(queryRange, codingRange)) {
               print(paste("OUTSIDE at ", rangeStr(queryRange), " -------------------", sep=''))
          } else {
               print(paste("Partially inside at ", rangeStr(queryRange), " ****", sep=''))
          }

          qCode    <- IRanges::intersect(queryRange, codingRange) # only the part of the query that is ALSO in the coding region
          uCode    <- IRanges::setdiff(queryRange, codingRange) # subtract out the coding part from the query

          xprint("Coding part of that: ", rangeStr(qCode))
          xprint("UTR part of that: ", rangeStr(uCode))
          
          #if (queryRange$START
     }
     
}

#for (i in seq_along(column9)) {

# Now let's go through them one TRANSCRIPT at a time











