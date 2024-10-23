#!/bin/bash
set -eu

# Script to summarise 5'P read ends (e.g. PARE-seq, GMUCT) across features of interest in bedGraph format and scale to reads per million (RPM)

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_bedgraph_5p.sh <.BAM> <layout: SE/PE> <bedfile annotation> <feature name> <distance>"
echo "EXAMPLE: BAM_to_bedgraph_5p.sh col0_rep1.sorted.bam PE Arabidopsis_thaliana.TAIR10.54_stop.bed stop 50"
exit 1
fi

smp=$1
lay=$2
bedfile=$3
out=$4
dis=$5

echo ""
echo "sample = $1"
echo "layout = $2"
echo "bedfile = $3"
echo "feature = $4"
echo "distance = $5"
echo ""

echo "calculate scaling factor"
if [[ "$lay" == "SE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)"); fi
if [[ "$lay" == "PE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)/2"); fi

echo "BAM to bed..."
bedtools genomecov -bg -5 -scale $scl -ibam $smp > ${smp%%bam}5p.bed

echo 'bedtools for coverage across chosen features...'
closestBed -D "b" -a ${smp%%bam}5p.bed -b $bedfile > ${smp%%.bam}_${out}.5p.bed

echo 'subset ...'
awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.5p.bed > ${smp%%.bam}_${out}_${dis}bp.5p.bed

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}_${dis}bp.5p.bed
	
echo 'cleaning'
rm -v ${smp%%bam}5p.bed ${smp%%.bam}_${out}.5p.bed ${smp%%.bam}_${out}_${dis}bp.5p.bed


