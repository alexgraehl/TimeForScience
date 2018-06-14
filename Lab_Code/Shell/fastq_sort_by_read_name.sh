#!/bin/bash
set -euo pipefail

INFQ=$1

OUTFQ=${INFQ}
OUTFQ=${OUTFQ/.fq/.sorted.fq}
OUTFQ=${OUTFQ/.fastq/.sorted.fastq}

if [[ "${INFQ}" == "${OUTFQ}" ]]; then
    echo "Something went wrong... the input and output filenames are the same???"
    OUTFQ="${OUTFQ}.sorted.fastq.gz"
fi
zcat "${INFQ}" | paste - - - - | sort -k1,1g -t " " | tr "\t" "\n" | gzip > "${OUTFQ}"

(>&2 echo "[STDERR] Generated the output file --> ${OUTFQ}")

