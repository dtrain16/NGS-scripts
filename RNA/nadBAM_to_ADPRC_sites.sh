#!/bin/bash
set -u

# Compute NAD enrichment per site (start of read) by providing +ADPRC and -ADPRC samples as BAM aligned reads

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots


if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: graft_nad_adprc_enrich.sh <+adprc bam> <-adprc bam> <sample name> <bedfile> <feature>"
echo "EXAMPLE: graft_nad_adprc_enrich.sh WT_plus-a_rep1.sorted.bam WT_minus-a_rep1.sorted.bam WT_rep1 tss.bed tss"
exit 1
fi

smp_p=$1
smp_m=$2
out=$3
bedfile=$4
feature=$5

echo "###############"
echo "+ADPRC sample = $1"
echo "-ADPRC sample= $2"
echo "Sample ID = $3"
echo "Bedfile = $4"
echo "Feature name = $5"
echo "#############"

echo "calculate scaling factors"
scl_p=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp_p)") 
scl_m=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp_m)")

echo "BAM to bed..."
bedtools genomecov -bg -5 -scale $scl_p -ibam $smp_p > ${smp_p%%.bam}.5p.bed
bedtools genomecov -bg -5 -scale $scl_m -ibam $smp_m > ${smp_m%%.bam}.5p.bed

echo "Combine + and -ADPRC samples to calculate per-nt NAD%"
bedtools unionbedg -names plus minus -i ${smp_p%%.bam}.5p.bed ${smp_m%%.bam}.5p.bed | awk 'BEGIN {FS=OFS="\t"} {prop = $4 / ($4 + $5) } {print $0, prop}' > ${out}.nad.5p.bed

echo "closestBed..."
closestBed -D "a" -a ${out}.nad.5p.bed -b $bedfile > ${out}_${feature}_nad.5p.bed
awk -F$'\t' '$NF<51 && $NF>-51' ${out}_${feature}_nad.5p.bed > ${out}_${feature}_20bp_nad.5p.bed 


echo "Maths ..."
Rscript /home/dganguly/scripts/RNA/rel_expression_plots_nad.r ${out}_${feature}_20bp_nad.5p.bed

echo 'cleaning'
rm -v ${smp_p%%.bam}.5p.bed ${smp_m%%.bam}.5p.bed ${out}.nad.5p.bed ${out}_${feature}_nad.5p.bed
