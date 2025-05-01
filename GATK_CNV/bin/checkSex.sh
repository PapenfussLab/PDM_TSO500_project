#!/bin/bash

BAM="$1"
sample="$2"
bed="$3"

module load samtools/1.19.2
samtools idxstats ${BAM} >${sample}.chromdepths

x_start=$(grep "chrX" ${bed} | cut -f2 | paste -sd+ | bc)
x_end=$(grep "chrX" ${bed} | cut -f3 | paste -sd+ | bc)
x_length=$((x_end - x_start))
x_map=`grep "chrX" ${sample}.chromdepths | awk '{print $3}'`
x_cov=$((100 * x_map / x_length ))

y_start=$(grep "chrY" ${bed} | cut -f2 | paste -sd+ | bc)
y_end=$(grep "chrY" ${bed} | cut -f3 | paste -sd+ | bc)
y_length=$((y_end - y_start))
y_map=`grep "chrY" ${sample}.chromdepths | awk '{print $3}'`
y_cov=$((100 * y_map / y_length))

xy_ratio=$(( x_cov / y_cov ))

if [[ ${xy_ratio} -gt 2 || "${y_cov}" == 0 ]]
then
sex="F"
else
sex="M"
fi

echo ${sex} > ${sample}_sex.txt

# echo "$sample"
# echo "X: ${x_cov}" >> ${sample}_sex.txt
# echo "Y: ${y_cov}" >> ${sample}_sex.txt
# echo "X/Y ratio: ${xy_ratio}" >> ${sample}_sex.txt
# echo "sex: ${sex}" >> ${sample}_sex.txt