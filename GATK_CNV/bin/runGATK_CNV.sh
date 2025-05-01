#!/bin/bash

sample="$1"
BAM_FILE="$2"
sex="$3"
INTERVAl="$4"
REF="${5}/genome.fa"

module load gatk/4.2.5.0

mkdir ${sample}_GATK

# Collect read counts
gatk CollectReadCounts \
-I ${BAM_FILE} \
-L ${INTERVAl} \
--interval-merging-rule OVERLAPPING_ONLY \
-O ${sample}_GATK/${sample}.counts.hdf5

# Denoise read counts
gatk --java-options "-Xmx8g" DenoiseReadCounts \
-I ${sample}_GATK/${sample}.counts.hdf5 \
--count-panel-of-normals read_count_PON_${sex}.hdf5 \
--standardized-copy-ratios ${sample}_GATK/${sample}_sdCR.tsv \
--denoised-copy-ratios ${sample}_GATK/${sample}_dnCR.tsv

## Collect allelic counts
gatk CollectAllelicCounts \
-I ${BAM_FILE} \
-R ${REF} \
-L ${INTERVAl} \
--interval-merging-rule OVERLAPPING_ONLY \
-O ${sample}_GATK/${sample}.allelicCounts.tsv

## Run Model segment
gatk --java-options "-Xmx8g" ModelSegments \
--denoised-copy-ratios ${sample}_GATK/${sample}_dnCR.tsv \
--allelic-counts ${sample}_GATK/${sample}.allelicCounts.tsv \
--output-prefix ${sample} \
-O ${sample}_GATK
