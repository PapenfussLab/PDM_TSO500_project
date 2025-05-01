#!/bin/bash
module load R/4.4.1
module load gatk/4.2.5.0

sex="$1"
export R_LIBS="$2" 
intervalPureCN="$3"
gvcfDir="$4"

PURECN="${R_LIBS}/PureCN/extdata"
mkdir normal_${sex}

if [[ "$sex" == "F" ]]; then
    gVCF="${gvcfDir}/TSO500_female.gvcf.gz"
elif [[ "$sex" == "M" ]]; then
    gVCF="${gvcfDir}/TSO500_male.gvcf.gz"
else
    echo "Error: unknown sex"
fi

ls *.bam | cat > normal_bam.list

Rscript $PURECN/Coverage.R \
    --out-dir normal_${sex} \
    --bam normal_bam.list \
    --intervals ${intervalPureCN} \
    --cores 8

ls normal_${sex}/*_loess.txt.gz | cat > normal_coverages.list

Rscript $PURECN/NormalDB.R \
    --out-dir normal_${sex} \
    --coverage-files normal_coverages.list \
    --normal-panel ${gVCF} \
    --genome hg19 \
    --assay TSO500