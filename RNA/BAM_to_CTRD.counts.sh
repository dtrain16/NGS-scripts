#!/bin/bash
set -u

# Determine counts within CTRD window (16-17 nt upstream from stop codons)
# Output based on featureCounts file format

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_CTRD.counts.sh <.BAM> <bedfile annotation>"
echo "EXAMPLE: BAM_to_CTRD.counts.sh WT_rep1.bam Arabidopsis_thaliana.TAIR10.54_stop.bed"
echo "annotation should provide coordinates of stop codons"
exit 1
fi

smp=$1
bedfile=$2

echo "sample = $1"
echo "bedfile = $2"

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_stop.5p.bed
awk -F$'\t' '$NF<-15 && $NF>-18' ${smp%%.bam}_stop.5p.bed | sort -k8,8 > ${smp%%.bam}_stop_CTRD.5p.bed 
groupBy -i ${smp%%.bam}_stop_CTRD.5p.bed -g 8,10 -c 4 -o sum > ${smp%%.bam}_CTRD.5p.counts

echo 'cleaning'
rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_stop.5p.bed ${smp%%.bam}_stop_CTRD.5p.bed

