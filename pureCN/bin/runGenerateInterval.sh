#!/bin/bash
export export R_LIBS="$1"
PURECN="${1}/PureCN/extdata"
REF="${2}/genome.fa"
MAPPABILITY="${3}"
BEDFILE="${4}"

module load R/4.4.1
module load gatk/4.2.5.0
module load samtools/1.21

Rscript $PURECN/IntervalFile.R \
    --in-file ${BEDFILE} --fasta ${REF} --out-file TSO500_intervals.txt \
    --genome hg19 \
    --mappability ${MAPPABILITY}