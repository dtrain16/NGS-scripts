#!/bin/bash
set -eu

# Multi-purpose read alignment pipeline with STAR
# Retrieve TAIR10 reference
# wget ftp://ftp.ensemblgenomes.org/pub/release-47/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# fasta index
# samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
# chromosome lengths
# cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len
# Build genome index
# STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbGTFfile Arabidopsis_thaliana.TAIR10.54.gtf --genomeDir /path/to/GenomeDir/ --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

### CONDA environment is installed
# conda create --name STAR_v1
# conda install -n STAR_v1 -c bioconda fastqc
# conda install -n STAR_v1 -c bioconda star
# conda install -n STAR_v1 -c grst trim_galore ## outdated version, install manually

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: STAR_pipe_v1.sh <SE/PE> <fastq R1> <R2> </path/to/index> <fileID>"
echo "EXAMPLE: STAR_pipe_v1.sh SE sample.fastq /home/dganguly/ref_seqs/STAR/ sample_rep1"
exit 1
fi

###
# SINGLE END
###

if [ "$1" == "SE" ]; then

# requirements
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: STAR_pipe_v1.sh <SE> <R1> </path/to/index> <fileID>"
echo "EXAMPLE: STAR_pipe_v1.sh SE sample.fastq /home/dganguly/ref_seqs/STAR/ sample_rep1"
exit 1
fi

#gather input variables
type=$1
fq=$2;
index=$3; #path to subread indexed reference genome
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing single-end alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_star_${dow}
mv $fq ${fileID}_star_${dow}
cd ${fileID}_star_${dow}

# initial fastqc
mkdir 1_fastqc
fastqc -t 8 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

echo "Read trimming... "

# read trimming with trimgalore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --fastqc --fastqc_args "-t 8" ../$fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv $fq 0_fastq

# read alignment
mkdir 3_align
mv 2_read_trimming/${fq%%.fastq*}_trimmed.fq* -t 3_align/
cd 3_align

echo "Beginning alignment ..."

STAR --runThreadN 8 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn ${fq%%.fastq*}_trimmed.fq* --outFileNamePrefix $fileID --outSAMtype BAM SortedByCoordinate | tee -a  ../${fileID}_${dow}.log

echo "cleaning..."

outbam="${fileID}*.sortedByCoord.out.bam"
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
mv *trimmed.fq.gz ../2_read_trimming/

echo "Alignment complete"

fi

####
# PAIRED END
####

if [ "$1" == "PE" ]; then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: STAR_pipe_v1.sh <PE> <R1> <R2> </path/to/index> <fileID>"
echo "EXAMPLE: STAR_pipe_v1.sh PE sample_R1.fq sample_R2.fq /home/dganguly/STAR sample_rep1"
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
echo "Performing paired-end alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq1 $fq2"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_star_${dow}
mv $fq1 ${fileID}_star_${dow}
mv $fq2 ${fileID}_star_${dow}
cd ${fileID}_star_${dow}

# initial fastqc
mkdir 1_fastqc
fastqc -t 8 $fq1 $fq2 2>&1 | tee -a ${fileID}_${dow}.log
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

echo ""
echo "Performing adapter and low-quality read trimming... "
echo ""

# adapter and quality trimming with trim_galore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --fastqc --fastqc_args "-t 8" --paired ../$fq1 ../$fq2 2>&1 | tee -a ../${fileID}_${dow}.log
cd ../

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# STAR align
mkdir 3_align
mv 2_read_trimming/${fq1%%.fastq*}_val_1.fq* -t 3_align/
mv 2_read_trimming//${fq2%%.fastq*}_val_2.fq* -t 3_align/
cd 3_align/

# subjunc read alignment
STAR --runThreadN 8 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn ${fq1%%.fastq*}_val_1.fq*, ${fq2%%.fastq*}_val_2.fq* --outFileNamePrefix $fileID --outSAMtype BAM SortedByCoordinate | tee -a  ../${fileID}_${dow}.log

echo "cleaning..."

mv *_val_1.fq.gz ../2_read_trimming/
mv *_val_2.fq.gz ../2_read_trimming/

echo "Alignment complete"

fi

