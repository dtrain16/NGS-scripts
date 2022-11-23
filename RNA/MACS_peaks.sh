#!/bin/bash
set -eu

## Detect peaks from m6A-RNA-seq and test for differential enrichment

### CONDA environment is installed
# conda create --name <name>
# conda install -n <name> -c bioconda bedtools
# conda install -n <name> -c bioconda macs2

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: MACS_peaks.sh <read layout> <stranded> <control> <IP>"
echo "EXAMPLE: MACS_peaks.sh SE/PE un/fr/rf sample1-input_merged sample1-m6_rep1"
exit 1
fi


#gather input variables
type=$1
input=$2
ip=$3

if [[ "$lay" == "SE" ]] && [[ "$str"  == "un" ]] ; then

input_bam=${input}*.sortedByCoord.out.bam
ip_bam=${ip}*.sortedByCoord.out.bam

macs2 callpeak --nomodel --extsize 50 -c $input -t $ip -f BAM -g 32542107 -n ${ip%%}10B_m6_1_C_0hr_MACS.plus -q 1e-2

fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "rf" ]] ; then

input_bam1=${input}*.sortedByCoord.out.minus.bam
input_bam2=${input}*.sortedByCoord.out.plus.bam

ip_bam1=${ip}*.sortedByCoord.out.minus.bam
ip_bam2=${ip}*.sortedByCoord.out.plus.bam

macs2 callpeak --nomodel --extsize 50 -c $input -t $ip -f BAM -g 32542107 -n ${ip%%}10B_m6_1_C_0hr_MACS.plus -q 1e-2

fi


