#!/bin/bash

sex="$1"
pondb="$2"
REF="${3}/genome.fa"

module load gatk/4.2.5.0

gatk SelectVariants \
-R ${REF} \
-V gendb://${pondb} \
-O ${sex}.gvcf.gz
