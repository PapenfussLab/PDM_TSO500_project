#!/bin/bash

sample="$1"
bam="$2"
interval="$3"
REF="${4}/genome.fa"
gnomAD_VCF="${5}/af-only-gnomad.raw.sites.vcf"

module load gatk/4.2.5.0

gatk Mutect2 -R ${REF} \
-I ${bam} \
-L ${interval} \
--max-mnp-distance 0 \
--germline-resource ${gnomAD_VCF} \
--genotype-germline-sites true \
--genotype-pon-sites true \
--native-pair-hmm-threads 8 \
-O ${sample}_normal_mutect2.vcf