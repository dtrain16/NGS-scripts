#!/bin/bash
set -eu

# StringTie2 pipeline: assemble transcripts using aligned reads (BAM output e.g. from Subjunc or Hisat2), merge assemblies into combined annotation, and quantify transcript abundance (TPM) across each sample in the dataset. 
# software required: samtools, stringtie, gffcompare.

## e.g. generate list of sample names 
## dir *bam > files.txt

if [ "$#" -lt 2 ]; then
echo "Missing required arguments!"
echo "USAGE: stringtie_pipe_v1.sh <sample list> <library strandedness> <reference GFF3/GTF>"
echo "EXAMPLE: stringtie_pipe_v1.sh files.txt <un / fr / rf> Arabidopsis_thaliana.TAIR10.46.gff3"
exit 1
fi

## gather input variables
smpls=$(cat $1)
type=$2
ref=$3
dow=$(date +"%F-%H-%m")

echo ""
echo "Organising BAM files"
echo ""

## ensure all BAM files are sorted and indexed
new_smpls=""

for i in $smpls; 
	do
	if [[ $i != *sorted.bam ]]; then 
		samtools sort -@ 4 $i -o "${i%%.bam}.sorted.bam"
		fq="${i%%.bam}.sorted.bam"
		samtools index $fq		
		new_smpls="$new_smpls $fq";
	else new_smpls="$smpls";
	fi;
done

echo ""
echo "Transcript assembly per sample"
echo ""

## Use each sample to produce transcript assembly

if [ "$type" == "un" ] && [ -z "$ref" ]; then

echo ""
echo "Unstranded library with no reference annotation"
echo ""

for i in $new_smpls; do
	stringtie $i -o "${i%%.sorted.bam}_stringtie.out";
	done
fi

if [ "$type" == "un" ] && [ ! -z "$ref" ]; then

echo ""
echo "Unstranded library with annotation $ref"
echo ""

for i in $new_smpls; do
        stringtie $i -G $ref -o "${i%%.sorted.bam}_stringtie.out";
	done
fi

if [ "$type" == "fr" ] && [ -z "$ref" ]; then

echo ""
echo "Forward stranded without reference annotation"
echo ""

for i in $new_smpls; do
	stringtie $i --fr -o "${i%%.sorted.bam}_stringtie.out";
	done
fi

if [ "$type" == "fr" ] && [ ! -z "$ref" ]; then

echo ""
echo "Forward stranded with reference $ref"
echo ""

for i in $new_smpls; do
	stringtie $i --fr -G $ref -o "${i%%.sorted.bam}_stringtie.out";
	done
fi

if [ "$type" == "rf" ]  && [ -z "$ref" ]; then

echo ""
echo "Reverse stranded with reference annotation"
echo ""

for i in $new_smpls; do
	stringtie $i --rf -o "${i%%.sorted.bam}_stringtie.out";
	done
fi

if [ "$type" == "rf" ]  && [ ! -z "$ref" ]; then

echo ""
echo "Reverse stranded with reference $ref"
echo ""

for i in $new_smpls; do
        stringtie $i --rf -G $ref -o "${i%%.sorted.bam}_stringtie.out";
        done
fi

## Merge assemblies
echo ""
echo "Merging assemblies"
echo ""

strng="*_stringtie.out"

if [ ! -z "$ref" ]; then 
	stringtie --merge -i $strng -G $ref -o "merged_stringtie_out.gtf"; 
	else stringtie --merge -i $strng -o "merged_stringtie_out.gtf";
fi

## clean-up
rm *_stringtie.out

## Determine abundance of assembled transcript in each sample

echo ""
echo "Abundance estimation"
echo ""

if [ "$type" == "un" ]; then

for i in $new_smpls; do
	stringtie $i -e -G merged_stringtie_out.gtf -A "${i%%.sorted.bam}_abund.tab" -o "${i%%.sorted.bam}_stringtie_out.gtf";
	done
fi

if [ "$type" == "fr" ]; then

for i in $new_smpls; do
        stringtie $i -e --fr -G merged_stringtie_out.gtf -A "${i%%.sorted.bam}_abund.tab" -o "${i%%.sorted.bam}_stringtie_out.gtf";
        done
fi

if [ "$type" == "rf" ]; then

for i in $new_smpls; do
        stringtie $i -e --rf -G merged_stringtie_out.gtf -A "${i%%.sorted.bam}_abund.tab" -o "${i%%.sorted.bam}_stringtie_out.gtf";
        done
fi

echo "cleaning ..."

for i in $smpls;
        do
        if [[ $i != *sorted.bam ]]; then
                rm "${i%%.bam}.sorted.bam"
                rm "${i%%.bam}.sorted.bam.bai"
        fi;
done

mkdir abundance_estimates
mv *abund.tab -t abundance_estimates/

if [ ! -z "$ref" ]; then

echo ""
echo "compare with $ref"
echo ""

gffcompare -R -r $ref -o strtcmp merged_stringtie_out.gtf

mkdir gffcompare_results
mv strtcmp* -t gffcompare_results/

fi

echo "Complete!"
##############################


