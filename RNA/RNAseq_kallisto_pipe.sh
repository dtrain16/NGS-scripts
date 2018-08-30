#!/bin/bash

# Use kallisto to perform k-mer based transcript quantification
# https://www.nature.com/articles/nbt.3519

set -eu

if [ "$#" -lt 5 ]; then
	echo "Missing arguments!"
	echo "USAGE: kallisto.sh <SE,PE> <R1> <R2> <strandedness> <annotation> <name>"
	echo "strand: unstranded, fr_stranded, rf_stranded"
	echo "EXAMPLE: kallisto.sh PE SRR5724597_1.fastq.gz SRR5724597_2.fastq.gz unstranded AtRTD2_19April2016.fa col0-r1"
exit 1
fi

dow=$(date +"%F-%H-%m-%S")

###########
### SINGLE END
###########

if [ "$1" == "SE" ]; then
	# requirements
	if [ "$#" -ne 5 ]; then
		echo "Missing required arguments for single-end!"
		echo "USAGE: kallisto.sh <SE> <R1> <strandedness> <annotation> <name>"
		exit 1
	fi

R1=$2
strand=$3
annotation=$4
name=$5

kallito index -i 

fi

##########################
############# PAIRED END
##########################

if [ "$1" == "PE" ]; then
	# requirements
	if [ "$#" -ne 6 ]; then
		echo "Missing required arguments for single-end!"
		echo "USAGE: kallisto.sh <PE> <R1> <R2> <strandedness> <annotation> <name>" 
		exit 1
	fi

	# gather input variables
	R1=$2
	R2=$3
	strand=$4
	annotation=$5
	name=$6

echo "##################"
echo "Performing paired-end kallisto RNA-seq alignment"
echo "Type: $1"
echo "Input Files: $R1 $R2"
echo "Annotation: $annotation"
echo "Sample: $name"
echo "Time of analysis: $dow"
echo "##################"

# file structure
mkdir ${name}_kallisto_${dow}
mv $R1 $R2 -t ${name}_kallisto_${dow}
cd ${name}_kallisto_${dow}
mkdir 0_fastq
mv $R1 $R2 -t 0_fastq/

mkdir 1_fastqc 
fastqc -t 4 $R1 $R2 2>&1 | tee -a ${sample}_logs_${dow}.log
mv ${R1%%.fastq*}_fastqc* 1_fastqc
mv ${R2%%.fastq*}_fastqc* 1_fastqc

### Read trimming
echo "Adapter and quality trimming"
mkdir 2_trimmed_fastq
scythe -a $HOME/scripts/TruSeq-adapters.fa -p 0.1 $R1 > ${R1%%.fastq*}_noadapt.fastq 2>&1 | tee -a ${name}_logs_${dow}.log

scythe -a $HOME/scripts/TruSeq-adapters.fa -p 0.1 $R2 > ${R2%%.fastq*}_noadapt.fastq 2>&1 | tee -a ${name}_logs_${dow}.log

sickle pe -f ${R1%%.fastq*}_noadapt.fastq -r ${R2%%.fastq*}_noadapt.fastq -o ${R1%%.fastq*}_trimmed.fastq -p ${R2%%.fastq*}_trimmed.fastq -t sanger -q 20 -l 20 2>&1 | tee -a ${name}_logs_${dow}.log

rm *noadapt.fastq
mv *trimmed.fastq -t 2_trimmed_fastq

## FastQC
mkdir 3_trimmed_fastqc

fastqc -t 4 2_trimmed_fastq/${R1%%.fastq*}_trimmed.fastq* 2_trimmed_fastq/${R2%%.fastq*}_trimmed.fastq* 2>&1 | tee -a ${name}_logs_${dow}.log

mv *trimmed_fastqc* 3_trimmed_fastqc

## Generate kallisto index
echo "kallisto"

kallisto index -i "${annotation%%.fa}.idx" $annotation 2>&1 | tee -a ${name}_logs_${dow}.log

kallisto quant -i ${annotation%%.fa}.idx -t 4 --bias 2_trimmed_fastq/${R1%%.fastq*}_trimmed.fastq* 2_trimmed_fastq/${R2%%.fastq*}_trimmed.fastq* -o ./ 2>&1 | tee -a ${name}_logs_${dow}.log

mv abundance.hd5 ${name}.hd5
mv abundance.tsv ${name}.tsv

echo "Complete"

fi

