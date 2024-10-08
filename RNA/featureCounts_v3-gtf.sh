#!/bin/bash

# Perform featureCounts on aligned BAM files to assign counts across specified features in a GTF file

set -eu

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts_v3-gtf.sh <filename> <SE/PE> <library strandedness> <bedfile> <gtf feature>"
echo "EXAMPLE: RNAseq_featureCounts_v3-gtf.sh col0-r1.bam PE 2 Arabidopsis.gtf gene"
echo "library strandedness: 0 = unstranded, 1 = stranded, 2 = reverse stranded"
exit 1
fi

sample=$1
layout=$2
strand=$3
bedfile=$4
feat=$5
outname="${sample%%.bam*}_${feat}.gtf.counts"

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3 where 0 = unstranded, 1 = stranded, 2 = reverse stranded"
echo "Counts at $5 in $bedfile"
echo ""

if [[ $layout == "SE" ]]; then 

featureCounts -F GTF -C -T 4 -M -t $feat -g gene_id -O -s $strand -a $bedfile -o $outname $sample
	
fi
	
if [[ $layout == "PE" ]]; then 

featureCounts -F GTF -p	-C -M -T 4 -t $feat -g gene_id -O -s $strand -a $bedfile -o $outname $sample

fi

# -C - Do not count read pairs matching different chromosomes
# -t - Specify feature
# -g - Specify attribute
# -O - Keep reads assigned to multiple features
# -p - paired-end reads, count fragments
# -M - Multi-mapping reads will also be counted.

echo "DONE"
