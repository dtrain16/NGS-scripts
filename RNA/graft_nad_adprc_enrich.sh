#!/bin/bash
set -u

# Compute NAD enrichment per site (start of read) by providing +ADPRC and -ADPRC samples as BAM aligned reads

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots


if [ "$#" -lt 6 ]; then
echo "Missing arguments!"
echo "USAGE: graft_nad_adprc_enrich.sh <+adprc bam> <-adprc bam> <sample name>"
echo "EXAMPLE: graft_nad_adprc_enrich.sh WT_plus-a_rep1.sorted.bam WT_minus-a_rep1.sorted.bam WT_rep1"
exit 1
fi

smp_p=$1
smp_m=$2
bedfile=$3
out=$4

echo "###############"
echo "sample + = $1"
echo "sample - = $2"
echo "bedfile = $3"
echo "sample = $4"
echo "#############"

echo "calculate scaling factors"
scl_p=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp_p)") 
scl_m=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp_m)")

echo "BAM to bed..."
bedtools genomecov -bg -5 -scale $scl_p -ibam $smp_p > ${smp_p%%.bam}.5p.bed
bedtools genomecov -bg -5 -scale $scl_m -ibam $smp_m > ${smp_m%%.bam}.5p.bed

echo "Combine + and - ADPRC samples to find enriched sites"
bedtools unionbedg -header -names plus minus -i ${smp_p%%.bam}.5p.bed ${smp_m%%.bam}.5p.bed > ${out}.5p.bed


