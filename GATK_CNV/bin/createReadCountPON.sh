#!/bin/bash

sex=$1
IFS=',' read -r -a readcounts <<< "$2"
condaDir="$3"

module load miniconda3
conda activate $condaDir

# Create PON only works with older version of GATK. Installed gatk 4.1.2.0 in conda
PATH=/vast/ai_projects/melanoma_depigmentation_biomarker_discovery/peinan/conda/gatk4/bin:$PATH 
pon_command="gatk --java-options "-Xmx8000m" CreateReadCountPanelOfNormals"

for rc in ${readcounts[@]}
do 
    pon_command+=" -I ${rc}"
done

pon_command+=" -O read_count_PON_${sex}.hdf5 --number-of-eigensamples 10"
eval $pon_command