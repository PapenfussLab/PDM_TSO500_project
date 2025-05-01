#!/bin/bash

module load gatk/4.2.5.0

sample="$1" 
bam="$2" 
sex="$3" 
refDir="$4" 
intervalGATK="$5" 
ponvcfDir="$6"
gnomAD="$7"
gnomAD_VCF="${gnomAD}/af-only-gnomad.raw.sites.vcf"

mkdir ${sample}_mutect2

if [[ "$sex" == "F" ]]; then
    pon_VCF="${ponvcfDir}/TSO500_female.vcf.gz"
elif [[ "$sex" == "M" ]]; then
    pon_VCF="${ponvcfDir}/TSO500_male.vcf.gz"
else
    echo "Error: unknown sex"
fi


gatk Mutect2 -R ${refDir}/genome.fa \
-I ${bam} \
-L ${intervalGATK} \
--germline-resource ${gnomAD_VCF} \
--panel-of-normals ${pon_VCF} \
--genotype-germline-sites true \
--genotype-pon-sites true \
-O ${sample}_mutect2/${sample}_mutect2.vcf
