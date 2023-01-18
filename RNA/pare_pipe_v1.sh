#!/bin/bash
set -eu

# Read alignment for PARE libraries with STAR
# SE only, trim reads to 20 bp, and remove length cutoff for trimming

# Retrieve TAIR10 reference and prepare STAR index
# wget ftp://ftp.ensemblgenomes.org/pub/release-47/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
# cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len

# Build STAR genome index
# STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbGTFfile Arabidopsis_thaliana.TAIR10.54.gtf --genomeDir /path/to/GenomeDir/ --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

### CONDA environment is installed
# conda create --name <name>
# conda install -n <name> -c bioconda fastqc
# conda install -n <name> -c bioconda star
# conda install -n <name> -c bioconda bedtools
# conda install -c bioconda fastx_toolkit

if [ "$#" -lt 3 ]; then
echo "Missing required arguments!"
echo "USAGE: pare_pipe_v1.sh <fastq R1> </path/to/index> <fileID>"
echo "EXAMPLE: pare_pipe_v1.sh sample.fastq /home/dganguly/ref_seqs/STAR/TAIR10/GenomeDir sample_rep1"
exit 1
fi

#gather input variables
fq=$1
index=$2; #path to STAR index
fileID=$3;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing single-end alignment with the following parameters:"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_pare_${dow}
mv $fq ${fileID}_pare_${dow}
cd ${fileID}_pare_${dow}

## trim FASTQ file to 20bp read lengths (MMeI cut length)
if [[ $fq != *.gz ]];then fastx_trimmer -z -l 20 -i $fq -o ${fq%%.fastq}.20bp.fastq.gz; fi
if [[ $fq == *.gz ]];then zcat $fq | fastx_trimmer -z -l 20 -o ${fq%%.fastq.gz}.20bp.fastq.gz; fi

mkdir 0_fastq
mv $fq 0_fastq
fq=${fq%%.fastq.gz}.20bp.fastq.gz

# initial fastqc
mkdir 1_fastqc
fastqc -t 8 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

echo "Read trimming... "
# read trimming with trimgalore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --length 10 --fastqc --fastqc_args "-t 8" ../$fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mv $fq 0_fastq

# read alignment
mkdir 3_align
mv 2_read_trimming/${fq%%.fastq*}_trimmed.fq* -t 3_align/

cd 3_align
echo "Beginning alignment ..."
input=${fq%%.fastq*}_trimmed.fq*
STAR --runThreadN 8 --outFilterMismatchNmax 0 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn $input --outFileNamePrefix $fileID --outSAMtype BAM SortedByCoordinate | tee -a  ../${fileID}_logs_${dow}.log

echo "cleaning..."

outbam="${fileID}*.sortedByCoord.out.bam"
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
mv *trimmed.fq.gz ../2_read_trimming/

echo "Alignment complete"


