#!/bin/bash

# Use kallisto to perform k-mer based transcript quantification
# https://www.nature.com/articles/nbt.3519
# Build annotation index kallisto index -i annotation.idx -k 21 annotation.fa

set -eu

if [ "$#" -lt 5 ]; then
	echo "Missing arguments!"
	echo "USAGE: kallisto.sh <SE,PE> <R1> <R2> <strandedness> <index> <name>"
	echo "strand: unstranded, fr_stranded, rf_stranded"
	echo "EXAMPLE: kallisto.sh PE SRR5724597_1.fastq.gz SRR5724597_2.fastq.gz unstranded AtRTD2_19April2016.idx col0-r1"
exit 1
fi

dow=$(date +"%F")

###########
### SINGLE END
###########

if [ "$1" == "SE" ]; then
	# requirements
	if [ "$#" -ne 5 ]; then
		echo "Missing required arguments for single-end!"
		echo "USAGE: kallisto.sh <SE> <R1> <strandedness> <index> <name>"
		exit 1
	fi

type=$1
R1=$2
strand=$3
annotation=$4
name=$5

echo "##################"
echo "Performing single-end alignments with kallisto"
echo "Type: $type"
echo "Input Files: $R1"
echo "Annotation: $annotation"
echo "Sample: $name"
echo "Time of analysis: $dow"
echo "##################"

# file structure
mkdir ${name}_kallisto_${dow}
mv $R1 -t ${name}_kallisto_${dow}
cd ${name}_kallisto_${dow}

mkdir 0_fastq
mv $R1 -t 0_fastq/

### Read trimming & FastQC
echo "Read trimming and FastQC"

mkdir 1_trimmed_fastq
cd 1_trimmed_fastq
trim_galore --fastqc --fastqc_args "--threads 4" ../0_fastq/$R1 | tee -a ../${name}_logs_${dow}.log
cd ../

mkdir 2_quant/
mv 1_trimmed_fastq/*fq.gz 2_quant/
cd 2_quant/

## Kallisto quantification in unstranded, fr-stranded, and rf-stranded modes
echo "                      "
echo "kallisto"
echo "                      "

if [ $strand == "unstranded" ]; then

	kallisto quant -i $annotation -t 4 --bias --single ${R1%%.fastq*}_trimmed.fq* -b 10 -l 300 -s 100 -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

elif [ $strand == "fr_stranded" ]; then
        kallisto quant -i $annotation --fr-stranded -t 4 --bias --single ${R1%%.fastq*}_trimmed.fq* -b 10 -l 300 -s 100 -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

else kallisto quant -i $annotation --rf-stranded -t 4 --bias --single ${R1%%.fastq*}_trimmed.fq* -b 10 -l 300 -s 100 -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

fi

fi

##########################
############# PAIRED END
##########################

if [ "$1" == "PE" ]; then
	# requirements
	if [ "$#" -ne 6 ]; then
		echo "Missing required arguments for paired-end!"
		echo "USAGE: kallisto.sh <PE> <R1> <R2> <strandedness> <annotation> <name>" 
		exit 1
	fi

# gather input variables
type=$1
R1=$2
R2=$3
strand=$4
annotation=$5
name=$6

echo "##################"
echo "Performing paired-end alignments with kallisto"
echo "Type: $type"
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

### Read trimming & FastQC
echo "Read trimming and FastQC"

mkdir 1_trimmed_fastq
cd 1_trimmed_fastq
trim_galore --fastqc --fastqc_args "--threads 4" --paired ../0_fastq/$R1 ../0_fastq/$R2 | tee -a ../${name}_logs_${dow}.log
cd ../

mkdir 2_quant/
mv 1_trimmed_fastq/*fq.gz 2_quant/
cd 2_quant/

### kallisto quantification in unstranded, fr-stranded or rf-stranded modes
echo "                      "
echo "kallisto"
echo "                      "

if [ $strand == "unstranded" ]; then

	kallisto quant -i $annotation -t 4 --bias -b 10 ${R1%%.fastq*}_val*fq.gz ${R2%%.fastq*}_val*fq.gz -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

elif [ $strand == "fr_stranded" ]; then
	kallisto quant -i $annotation --fr-stranded -t 4 --bias -b 10 ${R1%%.fastq*}_val*fq.gz ${R2%%.fastq*}_val*fq.gz -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

else kallisto quant -i $annotation --rf-stranded -t 4 --bias -b 10 ${R1%%.fastq*}_val*fq.gz ${R2%%.fastq*}_val*fq.gz -o ./ 2>&1 | tee -a ../${name}_logs_${dow}.log

fi

fi

## archive trimmed fastq file
mv *fq.gz ../1_trimmed_fastq/

### extract quant files for local import
cd ../../
mkdir quant_files
fls=*_kallisto_${dow}
for i in $fls; do echo $i; test=${i%%_kallisto*}; mkdir quant_files/$test; cp ${i}/2_quant/* -t quant_files/$test; done

echo "complete"
