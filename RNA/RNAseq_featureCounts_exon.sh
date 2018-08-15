#!/bin/bash

# Perform featureCounts on aligned BAM files to assign transcript counts PER EXON based on the AtRTD2 reference transcript dataset (or quasi - for alternative splicing)
# Additional info https://www.biostars.org/p/321379/
# https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf 
# Zhang R et al (2017). A high quality Arabidopsis transcriptome for accurate transcript-level analysis of alternative splicing. Nucleic Acids Res. 45: 5061–5073.

set -eu

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts.sh <filename> <SE/PE> <library strandedness> <rtd2 or padded>"
echo "EXAMPLE: RNAseq_featureCounts.sh col0-r1.bam PE 2 rtd2"
echo "library strandedness: 0 = unstranded, 1 = stranded, 2 = reverse stranded"
exit 1
fi

sample=$1
layout=$2
strand=$3
ref=$4

if [[ $ref == "rtd2" ]]; then
	bedfile="$HOME/AtRTD2/AtRTD2_19April2016.gtf"
	outname="${sample%%.bam*}_RTD2.counts"
elif [[ $ref == "padded" ]]; then
	bedfile="$HOME/AtRTD2/AtRTDv2_QUASI_19April2016.gtf"
	outname="${sample%%.bam*}_quasi.counts"
else
	echo " bad argument - pick 'rtd2' or 'padded' "
	exit 1
fi

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3 where 0 = unstranded, 1 = stranded, 2 = reverse stranded"
echo "$bedfile"
echo ""
echo "Exon feature counting - $layout $strand featureCounts on $bedfile in $sample ..."
echo ""

if [[ $layout == "SE" ]]; then 

featureCounts -F GTF -C -T 4 -f -t exon -g gene_id -O -s $strand -a $bedfile -o $outname $sample 2>&1 | tee -a ../*log
	
fi
	
if [[ $layout == "PE" ]]; then 

featureCounts -F GTF -p	-C -T 4 -f -t exon -g gene_id -O -s $strand -a $bedfile -o $outname $sample 2>&1 | tee -a ../*log

fi

# -C - Do not count read pairs matching different chromosomes
# -f - Perform read counting at feature level (e.g. exons vs genes)
# -t - Specify meta-feature
# -g - Specify attribute
# -O - Keep reads assigned to multiple meta-features

echo "DONE"
