#!/bin/bash -u
set -e
set -o pipefail

# INPUT (argument to this program): A file that should be a 'subread featureCounts' counts file
# OUTPUT (to STDOUT): The headers and numeric data from the file, but not the annotation columns.

if [[ "$#" == 0 ]] ; then echo "[ERROR] You need to pass EXACTLY ONE 'subread featureCounts' matrix file into this script!"; exit 1; fi

FILE=$1
head -n 5 $FILE | grep '^#' | grep 'featureCounts' > /dev/null
if [[ $? == 0 ]]; then
    IS_FEATURE_COUNT_FILE=1
    #echo "Is this a feature count file? Answer is YES: $IS_FEATURE_COUNT_FILE" >&2
else
    IS_FEATURE_COUNT_FILE=0
    # >&2  means "echo to stderr"
    echo "This was apparently not a 'featureCounts' file, so we are exiting now!" >&2
    echo "[FATAL ERROR] FAILURE -- the input file ($FILE) was seemingly not a Subread featureCounts file!"
    exit 1
fi

cat $FILE | grep -v '^#' | cut -f 1,7- 
# Outputs just the plain matrix part of the file
