#!/bin/bash
set -u
set -e

# By Alex Williams, Jan 1, 2015. Runs "qstat -f" and selects important things out of the output, colorizes the output, and formats it into columns.
# Better than qstat -s because it shows THE ENTIRE JOB NAME which is impossible to get in qstat -s

HEAD=`mktemp`
YYY=`mktemp`
ZZZ=`mktemp`
PRECOLOR=`mktemp`
echo -e 'Job ID\tJob name\tUser\tStatus\tMax time\tElapsed time\tStart time\tSeconds Remaining' > $HEAD


qstat -f | perl -e 'my $jn=undef; while(<>) { chomp; if(m/Job_Name/i){$jn=$_;next;}; if(defined($jn)){ if(m/.*=.*/){print(qq{$jn\n});$jn=undef;} else {s/^\s+//;$jn.=$_;next;}; }   print qq{$_\n}; }' > $YYY
# Q: what is going on here? Answer: 'jn' is the job name, and it can sometimes WRAP to the next line, for example:

# Job Id: 91554.rigel.gladstone.internal
#  Job_Name = ns_620_nicole_stone_mm9_rna_follows_ns353_and_ns560_mar_16__s22
#             _cuffdiff__all                       <--  we want THIS part to be added onto the "..._mar16__s22" on the line above.
# Job_Owner = alexgw@rigel.gladstone.internal

# ... so what we do is look for Job_Name, and then we basically just remove newlines and whitespace from everything UNTIL we see something with an equals sign, which means we got to the next element.
#     This is because this data is super unstructured and hard to parse.

cat $YYY | egrep -i '(Job Id|Job_Name|Job_Owner|job_state|list.walltime)' | perl -pe 's/\.gladstone\.internal//g' | perl -pe 's/(.*Job_Owner.*)@.*/\1/i' | perl -pe 's/.*=\s//' | perl -pe 's/\n/\t/g' | perl -pe 's/(Job Id:\s|$)/\n/ig'                                                > $ZZZ

# Note: job name can be split across MULTIPLE lines!


qstat -f | egrep -i '(Job Id|used.walltime|walltime.remaining|start_time)'                                       | perl -pe 's/.*=\s//' | perl -pe 's/\n/\t/g' | perl -pe 's/(Job Id:\s|$)/\n/ig' | cut -f 2- | paste -d '' $ZZZ - | cat $HEAD - > $PRECOLOR

FIN="\e[0m" # "finish" coloring
UCOL="\e[1;33;41m" # Colored text for this user's job. Anything that literally matches $USER
UHID="\e[8m" # HIDDEN text for the spacer that shows up when the job belongs to someone who is NOT the user (in the ugly perl thing below)
cat $PRECOLOR \
    | perl -p -e 's/\tC\t/\tCOMPLETE\t/g' \
    | perl -p -e 's/\tE\t/\tEXITING \t/g' \
    | perl -p -e 's/\tQ\t/\tQUEUED  \t/g' \
    | perl -p -e 's/\tR\t/\tRUNNING \t/g' \
    | perl -p -e 's/\tH\t/\tON HOLD \t/g' \
    | perl -p -e       's/(.*Job ID.*)/\e[31m\t$1\t\e[0m/g' \
    | perl -p -e 's/(.*\bCOMPLETE\b.*)/\e[34m\t$1\t\e[0m/g' \
    | perl -p -e  's/(.*\bEXITING\b.*)/\e[34m\t$1\t\e[0m/g' \
    | perl -p -e   's/(.*\bQUEUED\b.*)/\e[33m\t$1\t\e[0m/g' \
    | perl -p -e  's/(.*\bRUNNING\b.*)/\e[32m\t$1\t\e[0m/g' \
    | perl -p -e  's/(.*\bON HOLD\b.*)/\e[35m\t$1\t\e[0m/g' \
    | perl -p -e  's/(.*Not Running.*)/\e[33m\t$1\t\e[0m/g' \
    | perl -p -e  "if (m/(.*)($USER)(.*)/) { s/(.*)($USER)(.*)/${UCOL} \$2 ${FIN}\\t\$1\$2\$3/ } else { my \$REST = q{ } x ${#USER}; s/^/${UHID} \$REST ${FIN}\\t/ }" \
    | column -s $'\t' -t

# Note: the last perl line above (with "s/($USER)/......" is to COLOR the CURRENT USERNAME if it's found, so you can easily see your jobs.

#    | perl -p -e  "s/(.*)($USER)(.*)/\e[33m\$2\e[0m\$1\$2\$3/g" 

