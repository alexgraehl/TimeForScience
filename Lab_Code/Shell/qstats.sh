#!/bin/bash
set -u
set -e

# By Alex Williams, Jan 1, 2015. Runs "qstat -f" and selects important things out of the output, colorizes the output, and formats it into columns.
# Better than qstat -s because it shows THE ENTIRE JOB NAME which is impossible to get in qstat -s

HEAD=`mktemp`
ZZZ=`mktemp`
PRECOLOR=`mktemp`
echo -e 'Job ID\tJob name\tUser\tStatus\tMax time\tElapsed time\tStart time\tSeconds Remaining' > $HEAD
qstat -f | egrep -i '(Job Id|Job_Name|Job_Owner|job_state|list.walltime)' | perl -pe 's/\.gladstone\.internal//g' | perl -pe 's/(.*Job_Owner.*)@.*/\1/i' | perl -pe 's/.*=\s//' | perl -pe 's/\n/\t/g' | perl -pe 's/(Job Id:\s|$)/\n/ig'                                                > $ZZZ
qstat -f | egrep -i '(Job Id|used.walltime|walltime.remaining|start_time)'                                       | perl -pe 's/.*=\s//' | perl -pe 's/\n/\t/g' | perl -pe 's/(Job Id:\s|$)/\n/ig' | cut -f 2- | paste -d '' $ZZZ - | cat $HEAD - > $PRECOLOR

cat $PRECOLOR \
    | perl -p -e 's/\tC\t/\tCOMPLETE\t/g' \
    | perl -p -e 's/\tE\t/\tEXITING \t/g' \
    | perl -p -e 's/\tQ\t/\tQUEUED  \t/g' \
    | perl -p -e 's/\tR\t/\tRUNNING \t/g' \
    | perl -p -e 's/\tH\t/\tON HOLD \t/g' \
    | perl -p -e       's/(.*Job ID.*)/\033[31m\t$1\t\033[0m/g' \
    | perl -p -e 's/(.*\bCOMPLETE\b.*)/\033[34m\t$1\t\033[0m/g' \
    | perl -p -e  's/(.*\bEXITING\b.*)/\033[34m\t$1\t\033[0m/g' \
    | perl -p -e   's/(.*\bQUEUED\b.*)/\033[33m\t$1\t\033[0m/g' \
    | perl -p -e  's/(.*\bRUNNING\b.*)/\033[32m\t$1\t\033[0m/g' \
    | perl -p -e  's/(.*\bON HOLD\b.*)/\033[35m\t$1\t\033[0m/g' \
    | perl -p -e  's/(.*Not Running.*)/\033[33m\t$1\t\033[0m/g' \
    | perl -p -e            "s/($USER)/\033[7;49;93m\$1\033[0m/g" \
    | column -s $'\t' -t

# Note: the last perl line above (with "s/($USER)/......" is to COLOR the CURRENT USERNAME if it's found, so you can easily see your jobs.

