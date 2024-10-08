#!/bin/bash
set -u

# Calculate terminal stalling index based on 5'-P end counts at STOP codon

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_TSI.sh <.BAM> <layout:SE/PE> <bedfile annotation>"
echo "EXAMPLE: BAM_to_TSI.sh col0_rep1.sorted.bam PE Arabidopsis_thaliana.TAIR10.54_stop.bed"
echo "annotation provided should represent coordinates of stop codons"
exit 1
fi

smp=$1
lay=$2
bedfile=$3

echo "sample = $1"
echo "layout = $2"
echo "bedfile = $3"

echo "calculate scaling factor"
if [[ "$lay" == "SE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)"); fi
if [[ "$lay" == "PE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)/2"); fi

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_stop.5p.bed
awk -F$'\t' '$NF<1 && $NF>-101' ${smp%%.bam}_stop.5p.bed > ${smp%%.bam}_stop_TSI.5p.bed 

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/TSI_calculation.r ${smp%%.bam}_stop_TSI.5p.bed

echo 'cleaning'
rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_stop.5p.bed


