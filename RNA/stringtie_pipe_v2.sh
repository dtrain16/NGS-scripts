#!/bin/bash
set -eu

# StringTie2 v2: use stringtie to compute transcript abundance for pre-assembled transcripts and export TPM.

## e.g. generate list of sample names 
## dir *bam > files.txt

if [ "$#" -lt 3 ]; then
echo "Missing required arguments!"
echo "USAGE: stringtie_pipe_v2.sh <sample list> <library strandedness> <reference GFF3/GTF>"
echo "EXAMPLE: stringtie_pipe_v2.sh files.txt <un / fr / rf> merged"
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

## Transcript abundance

echo ""
echo "Abundance estimation (TPM)"
echo ""

if [ "$type" == "un" ]; then

for i in $new_smpls; do
	stringtie $i -e -G $ref -A "${i%%.sorted.bam}_gene.tab" -o "${i%%.sorted.bam}_stringtie.gtf";
	rm "${i%%.sorted.bam}_stringtie_out.gtf";
	mv t_data.ctab ${i%%.sorted.bam}_tdata.ctab
	done
fi

if [ "$type" == "fr" ]; then

for i in $new_smpls; do
        stringtie $i -e --fr -G $ref -A "${i%%.sorted.bam}_gene.tab" -o "${i%%.sorted.bam}_stringtie.gtf";
        rm "${i%%.sorted.bam}_stringtie_out.gtf";
	mv t_data.ctab ${i%%.sorted.bam}_tdata.ctab
	done
fi

if [ "$type" == "rf" ]; then

for i in $new_smpls; do
        stringtie $i -e --rf -G $ref -A "${i%%.sorted.bam}_gene.tab" -o "${i%%.sorted.bam}_stringtie.gtf";
	rm "${i%%.sorted.bam}_stringtie_out.gtf";
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

echo "Complete!"
##############################

