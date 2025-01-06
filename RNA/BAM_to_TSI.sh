#!/bin/bash
set -u

# Calculate terminal stalling index based on 5'-P end counts at STOP codon

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_TSI.sh <.BAM> <bedfile annotation>"
echo "EXAMPLE: BAM_to_TSI.sh col0_rep1.sorted.bam Arabidopsis_thaliana.TAIR10.54_stop.bed"
echo "annotation provided should represent coordinates of stop codons"
exit 1
fi

smp=$1
bedfile=$2

echo "sample = $1"
echo "bedfile = $2"

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_stop.5p.bed
awk -F$'\t' '$NF<1 && $NF>-51' ${smp%%.bam}_stop.5p.bed > ${smp%%.bam}_stop_TSI.5p.bed 

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/TSI_calculation.r ${smp%%.bam}_stop_TSI.5p.bed

echo 'cleaning'
rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_stop.5p.bed


