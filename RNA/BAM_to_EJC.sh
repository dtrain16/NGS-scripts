#!/bin/bash
set -u

# Script to extract genome coverage across features of interest
# optimised to caluclate 5'-P end frequency at exons from PARE or GMUCT
# only SE

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_EJC.sh <.BAM> <bedfile annotation> <feature name>"
echo "treat as unstranded only (degradome)"
echo "EXAMPLE: BAM_to_EJC.sh col0_rep1.sorted.bam Arabidopsis_thaliana.TAIR10.54_exon.bed exon"
exit 1
fi

smp=$1
bedfile=$2
out=$3

echo ""
echo "sample = $1"
echo "strand = unstranded"
echo "bedfile = $2"
echo "feature = $3"
echo ""

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_${out}.5p.bed
awk -F$'\t' '$NF<2 && $NF>-2' ${smp%%.bam}_${out}.5p.bed > ${smp%%.bam}_${out}_10bp.5p.bed 

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc.r ${smp%%.bam}_${out}_10bp.5p.bed

echo 'cleaning'
rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_${out}.5p.bed


