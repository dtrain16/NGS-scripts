#!/bin/bash
set -eu

# Script to identify sites with max 5'P read ends (e.g. PARE-seq, GMUCT) at features of interest in bedGraph format

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_bedgraph_max5p.sh <.BAM> <bedfile annotation> <feature name> <distance>"
echo "EXAMPLE: BAM_to_bedgraph_max5p.sh col0_rep1.sorted.bam Arabidopsis_thaliana.TAIR10.54_stop.bed stop 50"
exit 1
fi

smp=$1
bedfile=$2
out=$3
dis=$4

echo ""
echo "sample = $1"
echo "bedfile = $2"
echo "feature = $3"
echo "distance = $4"
echo ""

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%bam}5p.bed

echo 'bedtools for coverage across chosen features...'
closestBed -D "b" -a ${smp%%bam}5p.bed -b $bedfile > ${smp%%.bam}_${out}.5p.bed

echo 'subset ...'
awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.5p.bed > ${smp%%.bam}_${out}_${dis}bp.5p.bed

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}_${dis}bp.5p.bed
	
echo 'cleaning'
rm -v ${smp%%bam}5p.bed ${smp%%.bam}_${out}.5p.bed ${smp%%.bam}_${out}_${dis}bp.5p.bed


