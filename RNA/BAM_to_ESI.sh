#!/bin/bash
set -u

# Calculate EJC stalling index based on 5'-P end counts upstream of exon-exon junctions

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_ESI.sh <.BAM> <bedfile annotation>"
echo "EXAMPLE: BAM_to_ESI.sh col0_rep1.sorted.bam Arabidopsis_thaliana.TAIR10.54_exon-mRNA.bed"
echo "annotation provided should represent exons"
exit 1
fi

smp=$1
bedfile=$2

echo ""
echo "sample = $1"
echo "bedfile = $2"

echo "BAM to bed..."
bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_ejc.5p.bed
awk -F$'\t' '$NF<2 && $NF>-2' ${smp%%.bam}_ejc.5p.bed > ${smp%%.bam}_ejc_ESI.5p.bed 

echo 'do maths'
Rscript /home/dganguly/scripts/RNA/ESI_calculation.r ${smp%%.bam}_ejc_ESI.5p.bed

echo 'cleaning'
rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_ejc.5p.bed


