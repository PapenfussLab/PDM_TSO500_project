#!/bin/bash

sex="${1}"
INTERVAL=${2}
REF="${3}/genome.fa"


module load gatk/4.2.5.0

vcfFiles=`find . -name "*.vcf"`
GenomicsDB_command="gatk GenomicsDBImport -R ${REF} -L ${INTERVAL} --genomicsdb-workspace-path pondb_${sex}"

for vcf in ${vcfFiles}
do 
GenomicsDB_command+=" -V $vcf"
done

eval $GenomicsDB_command