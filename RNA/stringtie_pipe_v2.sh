#!/bin/bash
set -eu

# StringTie2 v2: use stringtie to compute transcript abundance for pre-assembled transcripts and export TPM.	
## e.g. generate list of sample names 	
## dir *bam > files.txt	

if [ "$#" -lt 3 ]; then	
echo "Missing required arguments!"	
echo "USAGE: stringtie_pipe_v2.sh <sample list> <strandedness> <reference GFF3/GTF>"	
echo "EXAMPLE: stringtie_pipe_v2.sh files.txt <un/fr/rf> merged"	
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
		samtools sort -@ 6 $i -o "${i%%.bam}.sorted.bam"	
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
	echo $i $type
	stringtie $i -e -G $ref -o "${i%%.sorted.bam}_stringtie.gtf";	
	done	
fi	

if [ "$type" == "fr" ]; then	

for i in $i; do
	echo $new_smpls $type
        stringtie $i -e --fr -G $ref -o "${i%%.sorted.bam}_stringtie.gtf";	
	done	
fi	

if [ "$type" == "rf" ]; then	

for i in $i; do
	echo $new_smpls $type	
        stringtie $i -e --rf -G $ref -o "${i%%.sorted.bam}_stringtie.gtf";	
	done	
fi	

echo "cleaning ..."	

for i in $smpls;	
        do	
        if [[ $i != *sorted.bam ]]; then	
                rm "${i%%.bam}.sorted.bam"	
        fi;	
done	

## export TPM

for i in $new_smpls; do
        Rscript $HOME/scripts/RNA/stringtie_extract_tpm.r "${i%%.sorted.bam}_stringtie.gtf"
done

mkdir tpm_abundance
mv *_stringtie.tpm tpm_abundance

echo "Complete!"	
##############################

