#!/bin/bash

sample="$1"
GATKDir="$2"
REF_dict="${3}/genome.dict"

module load gatk/4.2.5.0

# Plot segments need R 3.6. Installed older version of R in conda
PATH=/vast/ai_projects/melanoma_depigmentation_biomarker_discovery/peinan/conda/r3.6/bin:$PATH

mkdir tmp
gatk PlotModeledSegments \
--denoised-copy-ratios ${GATKDir}/${sample}_dnCR.tsv \
--segments ${GATKDir}/${sample}.modelFinal.seg \
--sequence-dictionary ${REF_dict} \
--output-prefix ${sample} \
--point-size-copy-ratio 3.0 \
-O .
