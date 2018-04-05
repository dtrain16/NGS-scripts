#!/bin/bash
set -eu

# gDNA whole-genome sequencing pipeline; qc trim (fastqc, trim_galore!), align (subread), and index raw reads

# Make sure genome index has been built using subread
# split all.fasta into chromosomes: samtools faidx genome.fasta chrX > chrX.fasta
# subread-buildindex -o TAIR10_subread_index TAIR10_Chr1.fasta TAIR10_Chr2.fasta TAIR10_Chr3.fasta TAIR10_Chr4.fasta TAIR10_Chr5.fasta TAIR10_ChrC.fasta TAIR10_ChrM.fasta

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: gdna_wgs_pipe.v1.sh <SE,PE> <fastq R1> <R2> <subread indexed genome> <fileID output>"
echo "EXAMPLE: gdna_wgs_pipe.v1.sh SE sample.fastq $HOME/TAIR10/chromosomes/TAIR10_subread_index sample-r1"
exit 1
fi

###
# SINGLE END
###

if [ "$1" == "SE" ]; then

# requirements
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: gdna_wgs_pipe.v1.sh <SE> <R1> <subread indexed ref genome> <fileID output>"
echo "EXAMPLE: gdna_wgs_pipe.v1.sh SE sample.fastq $HOME/TAIR10/chromosomes/TAIR10_subread_index sample-r1"
exit 1
fi

#gather input variables
type=$1
fq=$2;
index=$3; #path to subread indexed reference genome
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing single-end gDNA-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_gDNA_wgs_${dow}
mv $fq ${fileID}_gDNA_wgs_${dow}
cd ${fileID}_gDNA_wgs_${dow}

if [[ $fq != *.gz ]];then
gzip $fq
fq="${fq}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

echo "Performing quality-based read trimming... "

mkdir 2_trimgalore
cd 2_trimgalore
trim_galore ../$fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
cd ../

# repeat fastqc
echo "FASTQC ..."

mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq%%.fastq*}_trimmed.fastq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq 0_fastq

# subread align
mkdir 4_subread-align
mv 2_trim_galore/${fq%%.fastq*}_trimmed.fastq -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

# -t 0 = RNA-seq -t 1 = genomic DNA seq

subread-align -T 4 -t 1 -i $index -r ${fq%%.fastq*}_trimmed.fastq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

if [[ $fq%%.fastq}* != *".gz" ]]; then gzip ${fq%%.fastq*}_trimmed.fastq; fi

echo "cleaning..."

tmpbam="${fileID}.bam"
outbam="${fileID}.sorted.bam"
samtools sort -m 2G ${tmpbam} -o $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
rm -v ${tmpbam}
mv *trimmed.fastq.gz ../2_trimgalore/

echo "Alignment complete"

fi

#### 
# PAIRED END
####

if [ "$1" == "PE" ]; then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: gdna_wgs_pipe.v1.sh <PE> <R1> <R2> <subread indexed genome> <fileID output>"
echo "EXAMPLE: gdna_wgs_pipe.v1.sh PE sample_R1.fastq sample_R2.fastq $HOME/TAIR10/chromosomes/TAIR10_subread_index sample-r1"
exit 1
fi

#gather input variables
type=$1
fq1=$2;
fq2=$3;
index=$4; #path to genome index
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing paired-end RNA-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq1 $fq2"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_gDNA_wgs_${dow}
mv $fq1 ${fileID}_gDNA_wgs_${dow}
mv $fq2 ${fileID}_gDNA_wgs_${dow}
cd ${fileID}_gDNA_wgs_${dow}

if [[ $fq1 != *.gz ]];then
gzip $fq1
fq1="${fq1}.gz"
fi

if [[ $fq2 != *.gz ]];then
gzip $fq2
fq2="${fq2}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc $fq1 $fq2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

echo "Performing quality-based read trimming... "

mkdir 2_trimgalore
cd 2_trimgalore
trim_galore --paired ../$fq1 ../$fq2 2>&1 | tee -a ${fileID}_logs_${dow}.log
cd ../

# repeat fastqc
echo "FASTQC ..."

# fastqc again
mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq1%%.fastq*}_trimmed.fastq 2_trimgalore/${fq2%%.fastq*}_trimmed.fastq 2>&1 | tee -a ${fileID}_logs_${dow}.log

mv 2_trimgalore/${fq1%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc
mv 2_trimgalore/${fq2%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# subread align
mkdir 4_subread-align
mv 2_trimgalore/${fq1%%.fastq*}_trimmed.fastq -t 4_subread-align/
mv 2_trimgalore/${fq2%%.fastq*}_trimmed.fastq -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

# -t 0 = RNA-seq -t 1 = genomic DNA seq

subread-align -T 4 -t 1 -i ${index} -r ${fq1%%.fastq*}_trimmed.fastq -R ${fq2%%.fastq*}_trimmed.fastq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

echo "cleaning..."

if [[ $fq1%%.fastq}* != *".gz" ]]; then gzip ${fq1%%.fastq*}_trimmed.fastq; fi
if [[ $fq2%%.fastq}* != *".gz" ]]; then gzip ${fq2%%.fastq*}_trimmed.fastq; fi

mv *trimmed.fastq.gz ../2_trimgalore/

echo "Alignment complete"

fi

