#!/bin/bash

module load R/4.4.1
module load gatk/4.2.5.0

sample="$1"
sex="$2"
mutect2PATH="$3"
gatkPATH="$4"
export R_LIBS="$5"


if [[ "$sex" == "F" ]]; then
    normalDB=$(find normal_F/ -name "normalDB*.rds")
    mappingbias=$(find normal_F/ -name "mapping_bias*.rds")
elif [[ "$sex" == "M" ]]; then
    normalDB=$(find normal_M/ -name "normalDB*.rds")
    mappingbias=$(find normal_M/ -name "mapping_bias*.rds")
else
    echo "Error: unknown sex"
fi

PURECN="${R_LIBS}/PureCN/extdata"

sample_short=$(echo "$sample" | sed 's/^PDM_//; s/^SCM_//')
hdf5_FILE=${gatkPATH}/${sample_short}.counts.hdf5
logR_FILE=${gatkPATH}/${sample_short}_dnCR.tsv 
seg_FILE=${gatkPATH}/segments/${sample}.modelFinal.seg
#vcf_FILE=${gatkPATH}/${sample_short}.allelicCounts.tsv

mkdir ${sample}_GATK_PureCN

if [[ "$sample_short" == "BX1701760A6" ]]; then
    echo "sample excluded"
else
    Rscript $PURECN/PureCN.R \
        --sampleid ${sample} \
        --out ${sample}_GATK_PureCN \
        --tumor ${hdf5_FILE} \
        --log-ratio-file ${logR_FILE} \
        --seg-file ${seg_FILE} \
        --mapping-bias-file ${mappingbias} \
        --vcf ${mutect2PATH}/${sample}_mutect2.vcf \
        --fun-segmentation Hclust \
        --max-non-clonal 0.5 \
        --max-ploidy 4 \
        --genome hg19 \ 
        --cores 8 \
        --force --post-optimize --seed 123
fi