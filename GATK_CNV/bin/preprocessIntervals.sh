#!/bin/bash

REF="${1}/genome.fa"
BED="$2"

module load gatk/4.2.5.0

gatk PreprocessIntervals \
-R ${REF} \
-L ${BED} \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--padding 100 \
-O preprocessed_hg19.interval_list
