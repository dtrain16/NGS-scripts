#!/bin/bash

# Use featureCounts to assign counts to annotated features (e.g. genes, transposons) from aligned BAM files

set -eu

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts.sh <filename> <SE/PE> <library strandedness> <bedfile> <bedfile format> <outname>"
echo "EXAMPLE: RNAseq_featureCounts.sh col0-r1.sorted.bam PE 1 AtRTD2_19April2016.gtf gtf RTD2"
echo "library strandedness: 0 = unstranded, 1 = stranded, 2 = reverse stranded"
echo "format: saf, bed, gtf"
exit 1
fi

sample=$1
layout=$2
strand=$3
bedfile=$4
format=$5
outname=$6

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "bedfile = $4 ($5) ($6)"
echo ""
echo "$layout $strand featureCounts on $bedfile $format ($outname) in $sample ..."
echo ""

if [[ $format == "saf" ]]; then
	
	if [[ $layout == "SE" ]]; then
		featureCounts\
			-F 'SAF'\
			-C\
			-T 2\
			-s $strand\
			-a temp2.saf\
		        -o "${1%%.bam*}_${outname}.counts"\
		        $sample 2>&1 | tee -a ../*log
	fi

	if [[ $layout == "PE" ]]; then
		featureCounts\
			-F SAF\
			-p\
			-C\
			-T 2\
			-s $strand\
			-a temp2.saf\
			-o "${1%%.bam*}_${outname}.counts"\
			$sample 2>&1 | tee -a ../*log
	fi

fi

if [[ $format == "bed" ]]; then
	## convert BED to SAF format
	awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $bedfile > temp.saf
	awk 'BEGIN {print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print}' temp.saf > temp2.saf
	
	if [[ $layout == "SE" ]]; then 
		featureCounts\
			-F SAF\
			-C\
			-T 2\
			-s $strand\
			-a temp2.saf\
			-o "${1%%.bam*}_${outname}.counts"\
			$sample 2>&1 | tee -a ../*log
	fi
		
	if [[ $layout == "PE" ]]; then
	       featureCounts\
		       -F SAF\
			-p\
			-C\
			-T 2\
			-s $strand\
			-a temp2.saf\
			-o "${1%%.bam*}_${outname}.counts"\
			$sample 2>&1 | tee -a ../*log
	fi

fi

if [[ $format == "gtf" ]]; then
	if [[ $layout == "SE" ]]; then
		featureCounts\
			-F GTF\
			-C\
			-T 4\
			-s $strand\
			-a $bedfile\
			-o "${1%%.bam*}_${outname}.counts"\
			$sample 2>&1 | tee -a ../*log
	fi
	
	if [[ $layout == "PE" ]]; then
		featureCounts\
			-F GTF\
			-p\
			-C\
			-T 2\
			-s $strand\
			-a $bedfile\
			-o "${1%%.bam*}_${outname}.counts"\
		        $sample 2>&1 | tee -a ../*log
	fi
fi

rm temp*.saf -v

echo "DONE"
