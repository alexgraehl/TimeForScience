#!/bin/bash
set -euo pipefail

INFQ=$1

OUTFQ=${INFQ}
OUTFQ=${OUTFQ/.fq/.sorted.fq/}
OUTFQ=${OUTFQ/.fastq/.sorted.fastq/}

zcat "${INFQ}" | paste - - - - | sort -k1,1g -t \\\" \\\" | tr \\\"\t\\\" \\\"\n\\\" | gzip > "${OUTFQ}"
