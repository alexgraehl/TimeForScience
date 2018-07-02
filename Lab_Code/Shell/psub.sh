# Initial version from Tev Dincer: https://gist.github.com/udincer/345d4767eda0e2c21849e6746953402b

if [[ "$#" != 1 ]] ; then echo -e "[:ERR:]] Please specify a command filename (e.g. 'yourcommand.bsub.sh')"; exit 99; fi

TASKS_FILE=$1
#         word count ......      ...remove leading whitespace from 'wc' output
N_TASKS=$(wc -l < ${TASKS_FILE} | awk '{ sub(/^[ \t]+/, ""); print }')
PRIORITY="medium" # default queue
LOGDIR=logs/$(basename ${TASKS_FILE})
mkdir -p ${LOGDIR}
echo "[:OK :] Created the output log directory '${LOGDIR}'"

echo "[:FYI:] Input: a single file with one command per line: $1"
if [[ ! -d "${LOGDIR}" ]] ; then echo -e "[:ERR:]] The log directory must already exist. It is expected to be here: ${LOGDIR}"; exit 99; fi

echo "[:FYI:] Each line is sent off as a (separate) paralellized array job."
echo "[:OK :] Submitting ${N_TASKS} tasks to queue '${PRIORITY}' from file '${TASKS_FILE}'."

cat <<CMD
#!/bin/bash
set -uo pipefail
#BSUB -q ${PRIORITY}
#BSUB -J "$(basename ${TASKS_FILE})_$(date +%s)_[1-${N_TASKS}]"
#BSUB -oo "logs/$(basename ${TASKS_FILE})/$(date +%s)_%J-%I.bsub.log.out"
#BSUB -eo "logs/$(basename ${TASKS_FILE})/$(date +%s)_%J-%I.bsub.log.err"
#BSUB -R "rusage[mem=1]" # Request 1 GB of RAM
#BSUB -n 6,12            # Request between 6 and 12 cores, preferably 12?
#BSUB -L /bin/bash

cmd_multi=\$(awk "NR==\$LSB_JOBINDEX" $TASKS_FILE)
eval \${cmd_multi}

CMD

