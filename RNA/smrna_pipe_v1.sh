#!/bin/bash
set -eu

# Read alignment for smRNA-seq libraries with STAR
# SE only, trim adapters and low quality bases (remove length cutoff for trimming), trim reads to first 25 bp. Map with STAR with 0 mismatches and min mapped length 17.

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
echo "USAGE: smrna_pipe_v1.sh <fastq R1> </path/to/index> <fileID>"
echo "EXAMPLE: smrna_pipe_v1.sh sample.fastq /home/dganguly/ref_seqs/STAR/TAIR10/GenomeDir sample_rep1"
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
mkdir ${fileID}_srna_${dow}
mv $fq ${fileID}_srna_${dow}
cd ${fileID}_srna_${dow}

# gzip if unzipped input file
if [[ $fq != *.gz ]];then gzip $fq; fq="${fq}.gz"; fi

# initial fastqc
mkdir 1_fastqc
fastqc -t 8 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

echo "Read trimming... "
# Trim_galore: remove adapters and low quality base-calls, retain reads as small as 15 bp, generate fastqc report on trimmed reads.
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --length 15 --fastqc --fastqc_args "-t 8" ../$fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv $fq 0_fastq

# STAR alignment: trim reads to 25 bp max, align to STAR with 0 mismatches and minimum mapped length 17 bp.
mkdir 3_align
mv 2_read_trimming/${fq%%.fastq*}_trimmed.fq.gz -t 3_align/
cd 3_align
echo "Beginning alignment ..."

zcat ${fq%%.fastq*}_trimmed.fq.gz | fastx_trimmer -z -l 25 -o ${fq%%.fastq}.25bp.trimmed.fq.gz
input=${fq%%.fastq}.25bp.trimmed.fq.gz

STAR --runThreadN 8 --outFilterMismatchNmax 0 --outFilterMatchNmin 17 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn $input --outFileNamePrefix $fileID --outSAMtype BAM SortedByCoordinate | tee -a  ../${fileID}_logs_${dow}.log

echo "cleaning..."

outbam="${fileID}*.sortedByCoord.out.bam"
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
mv *trimmed.fq.gz ../2_read_trimming/

echo "Alignment complete"


