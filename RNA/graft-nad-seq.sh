#!/bin/bash
set -eu

# Single-end read alignment for GRAFT-NAD-seq libraries with stringent STAR parameters

# Prepare STAR index based on TAIR10 reference
# wget ftp://ftp.ensemblgenomes.org/pub/release-47/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
# cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len

# Build STAR genome index
# STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbGTFfile Arabidopsis_thaliana.TAIR10.54.gtf --genomeDir /path/to/GenomeDir/ --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

### CONDA environment is installed
# conda create --name <name>
# conda install -n <name> -c bioconda fastqc
# conda install -n <name> -c bioconda cutadapt
# conda install -n <name> -c bioconda star

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: graft-nad-seq.sh <fastq R1> </path/to/index> <adapter_file> <fileID>"
echo "EXAMPLE: graft-nad-seq.sh sample.fastq ~/ref_seqs/STAR/TAIR10/GenomeDir ~/GRAFT-NAD-seq_adapters.fa sample_rep1"
exit 1
fi


#gather input variables
fq=$1
index=$2;
adap=$3; #path to subread indexed reference genome
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing single-end alignment with the following parameters:"
echo "Input Files: $fq"
echo "genome index: $index"
echo "adapter sequences: $adap"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_graft_${dow}
mv $fq ${fileID}_graft_${dow}
cd ${fileID}_graft_${dow}

# initial fastqc
mkdir 1_fastqc
fastqc -t 8 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log

if [[ $fq == *"fq.gz" ]]; then mv ${fq%%.fq*}_fastqc* 1_fastqc; else mv ${fq%%.fastq*}_fastqc* 1_fastqc; fi

echo "Read trimming... "

# read trimming with trimgalore
mkdir 2_read_trimming
cd 2_read_trimming
trim_galore --length 30 --cores 8 --adapter_fasta $adap --fastqc --fastqc_args "-t 8" ../$fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv $fq 0_fastq

# read alignment
mkdir 3_align

if [[ $fq == *"fq.gz" ]]; then mv 2_read_trimming/${fq%%.fq*}_trimmed.fq* -t 3_align/; else
        mv 2_read_trimming/${fq%%.fastq*}_trimmed.fq* -t 3_align/;fi

cd 3_align

echo "Beginning alignment ..."

if [[ $fq == *"fq.gz" ]]; then input=${fq%%.fq*}_trimmed.fq*; else input=${fq%%.fastq*}_trimmed.fq*; fi

STAR --runThreadN 8 --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn $input --outFileNamePrefix "${fileID}_" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 8000000000 | tee -a  ../${fileID}_logs_${dow}.log

echo "cleaning..."

outbam="${fileID}*.sortedByCoord.out.bam"
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
mv *trimmed.fq.gz ../2_read_trimming/

echo "Alignment complete"


