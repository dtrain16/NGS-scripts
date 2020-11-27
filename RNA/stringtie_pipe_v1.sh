#!/bin/bash
set -eu

# StringTie2 v1: assemble transcripts using aligned reads (BAM output), merge assemblies, and quantify transcript abundance (-eB) for each sample. 
# software required: samtools, stringtie, gffcompare.

## e.g. generate list of sample names 
## dir *bam > files.txt

if [ "$#" -lt 3 ]; then
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
		new_smpls="$new_smpls $fq";
	else new_smpls="$new_smpls $i";
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

## Merge assemblies with abundance filters
echo ""
echo "Merging assemblies"
echo ""

strng="*_stringtie.out"
tpm_filter=$(ls -1 $strng | wc -l)
tpm_filter=$(expr $tpm / 2)

if [ ! -z "$ref" ]; then 
	stringtie --merge -i $strng -G $ref -f 0.05 -T $tpm_filter $ -o "merged_stringtie_out.gtf"; 
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
	stringtie $i -eB -G merged_stringtie_out.gtf -o "${i%%.sorted.bam}_stringtie_out.gtf";
	mv t_data.ctab ${i%%.sorted.bam}_tdata.ctab
	done
fi

if [ "$type" == "fr" ]; then

for i in $new_smpls; do
        stringtie $i -eB --fr -G merged_stringtie_out.gtf -o "${i%%.sorted.bam}_stringtie_out.gtf";
        mv t_data.ctab ${i%%.sorted.bam}_tdata.ctab
	done
fi

if [ "$type" == "rf" ]; then

for i in $new_smpls; do
        stringtie $i -eB --rf -G merged_stringtie_out.gtf -o "${i%%.sorted.bam}_stringtie_out.gtf";
        mv t_data.ctab ${i%%.sorted.bam}_tdata.ctab
	done
fi

echo "cleaning ..."

for i in $smpls;
        do
        if [[ $i != *sorted.bam ]]; then
                rm "${i%%.bam}.sorted.bam"
        fi;
done

## coverage estimates
mkdir abundance_estimates
mv *tdata.ctab -t abundance_estimates/
rm *ctab

## export TPM

for i in $new_smpls; do
	Rscript $HOME/scripts/RNA/stringtie_extract_tpm.r "${i%%.sorted.bam}_stringtie_out.gtf"
	mv *tpm -t abundance_estimates/
done

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


