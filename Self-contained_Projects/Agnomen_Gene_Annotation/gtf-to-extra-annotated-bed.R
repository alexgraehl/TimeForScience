# options(error=recover)

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
resMat <- matrix(data="", nrow=nrow(tomato), ncol=NUM_INTERESTING_ITEMS_IN_RESULT_VEC)

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

# Now let's go through them one TRANSCRIPT at a time
