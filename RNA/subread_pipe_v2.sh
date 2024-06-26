#!/bin/bash
set -eu

# RNA-sequencing pipeline: align trimmed reads using Subjunc
# Retrieve TAIR10 reference
# wget ftp://ftp.ensemblgenomes.org/pub/release-47/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# fasta index
# samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
# chromosome lengths
# cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len
# Build subread index
# subread-buildindex -o TAIR10.47_subread_index Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: RNAseq_subread_v2.sh <SE/PE> <fastq R1> <R2> <subread index> <fileID>"
echo "EXAMPLE: RNAseq_subread_v2.sh SE sample.fastq TAIR10_subread_index sample_rep1"
exit 1
fi

###
# SINGLE END
###

if [ "$1" == "SE" ]; then

# requirements
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: RNAseq_subread_v2.sh <SE> <R1> <subread index> <fileID>"
echo "EXAMPLE: RNAseq_subread_v2.sh SE sample.fastq TAIR10_subread_index sample_rep1"
exit 1
fi

#gather input variables
type=$1
fq=$2;
index=$3; #path to subread indexed reference genome
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing single-end RNA-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_subread_${dow}
mv $fq ${fileID}_subread_${dow}
cd ${fileID}_subread_${dow}

if [[ $fq != *.gz ]];then
gzip $fq
fq="${fq}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc -t 4 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

echo "Performing adapter and low-quality read trimming... "

# read trimming with trimgalore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --fastqc --fastqc_args "-t 4" ../$fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv $fq 0_fastq

# subread align
mkdir 3_subjunc
mv 2_read_trimming/${fq%%.fastq*}_trimmed.fq* -t 3_subjunc/
cd 3_subjunc

echo "Beginning alignment ..."

# subjunc aligner
subjunc -T 6 -i $index -r ${fq%%.fastq*}_trimmed.fq* -o  "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

if [[ ${fq%%.fastq*}* != *.gz ]]; then gzip ${fq%%.fastq*}* ; fi

echo "cleaning..."

tmpbam="${fileID}.bam"
outbam="${fileID}.sorted.bam"
samtools sort -m 2G ${tmpbam} -o $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
rm -v ${tmpbam}
mv *trimmed.fq.gz ../2_read_trimming/

echo "Alignment complete"

fi

####
# PAIRED END
####

if [ "$1" == "PE" ]; then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: RNAseq_subread_v2.sh <PE> <R1> <R2> <subread index> <fileID>"
echo "EXAMPLE: RNAseq_subread_v2.sh PE sample_R1.fq sample_R2.fq TAIR10_subread_index sample_rep1"
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
mkdir ${fileID}_subread_${dow}
mv $fq1 ${fileID}_subread_${dow}
mv $fq2 ${fileID}_subread_${dow}
cd ${fileID}_subread_${dow}

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
fastqc -t 4 $fq1 $fq2 2>&1 | tee -a ${fileID}_${dow}.log
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

echo ""
echo "Performing adapter and low-quality read trimming... "
echo ""

# adapter and quality trimming with trim_galore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --fastqc --fastqc_args "-t 4" --paired ../$fq1 ../$fq2 2>&1 | tee -a ../${fileID}_${dow}.log
cd ../

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# subread align
mkdir 3_subjunc
mv 2_read_trimming/${fq1%%.fastq*}_val_1.fq* -t 3_subjunc/
mv 2_read_trimming//${fq2%%.fastq*}_val_2.fq* -t 3_subjunc/
cd 3_subjunc/

echo "Beginning alignment ..."

# subjunc read alignment
subjunc -T 6 -i $index -r ${fq1%%.fastq*}_val_1.fq* -R ${fq2%%.fastq*}_val_2.fq* -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_${dow}.log

echo "cleaning..."

if [[ ${fq1%%.fastq*}* != *.gz ]]; then gzip ${fq1%%.fastq*}* ; fi
if [[ ${fq2%%.fastq*}* != *.gz ]]; then gzip ${fq2%%.fastq*}* ; fi

mv *_val_1.fq.gz ../2_read_trimming/
mv *_val_2.fq.gz ../2_read_trimming/

echo "Alignment complete"

fi

