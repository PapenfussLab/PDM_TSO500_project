#!/bin/bash

sample="$1"
bam="$2"
interval="$3"

module load gatk/4.2.5.0

gatk CollectReadCounts \
-I $bam \
-L $interval \
--interval-merging-rule OVERLAPPING_ONLY \
--format TSV \
-O ${sample}_normal_counts.tsv
