#!/bin/bash -u
set -e
set -o pipefail

# INPUT (argument to this program): A file that should be a 'subread featureCounts' counts file
# OUTPUT (to STDOUT): The headers and numeric data from the file, but not the annotation columns.

# Note: if you pass in "--base" as the second argument, it removes all directory names from the header and ONLY keeps the basenames!

if [[ "$#" == 0 ]] ; then echo "[ERROR] You need to pass EXACTLY ONE 'subread featureCounts' matrix file, and one optional second argument '--base' into this script!"; exit 1; fi

FILE=$1

if [[ -z "${2:-}" ]]; then
    # variable is not set or is empty
    ONLY_BASENAMES=0
elif [[ "$2" == "--base" ]]; then
    ONLY_BASENAMES=1
else
    echo "The second argument must either be ABSENT (for 'do not do anything to the basenames' or --base (for 'clip input names and only show the 'basenames'')"
    exit 1
fi

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

TEMP1=`mktemp`
TEMPHEADER=`mktemp`

cat $FILE | grep -v '^#' | cut -f 1,7- >| $TEMP1

# Remove file paths from filenames!

#echo $ONLY_BASENAMES

# | perl -pe 's/[._][sb]am\b//g'

if [[ "$ONLY_BASENAMES" == "1" ]]; then
    head -n 1 $TEMP1 | perl -e 'use File::Basename; while (<>) { chomp; my @a = split(/\t/, $_); my @z = map { File::Basename::basename($_) } @a; print join("\t", @z) . qq{\n}; }' > $TEMPHEADER
else
    head -n 1 $TEMP1  > $TEMPHEADER
fi
tail -n +2 $TEMP1 | cat $TEMPHEADER -
# This goes to STDOUT.

#echo "wrote $TEMP1 and $TEMP2"
/bin/rm $TEMP1 $TEMPHEADER
# Outputs just the plain matrix part of the file
