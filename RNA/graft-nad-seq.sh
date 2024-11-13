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
# conda install -n <name> -c bioconda seqkit 

if [ "$#" -lt 3 ] || [ "$#" -gt 3 ]; then
echo "Missing required arguments!"
echo "USAGE: graft-nad-seq.sh <fastq R2> </path/to/index> <fileID>"
echo "EXAMPLE: graft-nad-seq.sh sample.fastq ~/ref_seqs/STAR/TAIR10/GenomeDir sample_rep1"
exit 1
fi

#gather input variables
fq=$1
index=$2;
fileID=$3;
dow=$(date +"%F-%H-%m")

echo "##################"
echo "Performing single-end alignment with the following parameters:"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_graft-nad_${dow}
mv $fq ${fileID}_graft-nad_${dow}
cd ${fileID}_graft-nad_${dow}

# initial fastqc
mkdir 1_fastqc

fastqc -t 12 $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
if [[ $fq == *"fq.gz" ]]; then mv ${fq%%.fq*}_fastqc* 1_fastqc; else mv ${fq%%.fastq*}_fastqc* 1_fastqc; fi

echo "Extract branch sequence and trim adapters"
mkdir 2_read_trimming
mv $fq 2_read_trimming
cd 2_read_trimming

## extract reads beginning with branch sequence (GCTTGTTGTG) with flexibility at first and last base
if [[ $fq == *"fq.gz" ]]; 
	then seqkit grep -j 12 -s -r -p "^.CTTGTTGT" $fq -o ${fq%%.fq*}_branch.fq.gz ;
	else seqkit grep -j 12 -s -r -p "^.CTTGTTGT" $fq -o ${fq%%.fastq*}_branch.fq.gz; 
fi

if [[ $fq == *"fastq.gz" ]]; then fq_branch=${fq%%.fastq*}_branch.fq.gz; else fq_branch=${fq%%.fq*}_branch.fq.gz; fi

## cutadapt to remove branch seq at 5' end of read and universal PCR primer at 3' end of read
cutadapt -g "^NCTTGTTGTB" -g "^NCTTGTTGTBB" -g "^NCTTGTTGTBBB" -g  "^NCTTGTTGTBBBG" -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" -e 0.2 -m 25 -o "${fq_branch%%.fq*}_trimmed.fq.gz" ${fq_branch} 2>&1 | tee -a ../${fileID}_logs_${dow}.log

fastqc -t 12 ${fq_branch%%.fq*}_trimmed.fq.gz 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv 2_read_trimming/$fq  0_fastq/
mv 2_read_trimming/$fq_branch 0_fastq/

# read alignment
mkdir 3_align
mv 2_read_trimming/${fq_branch%%.fq*}_trimmed.fq.gz -t 3_align/

cd 3_align

echo "Beginning alignment ..."
input=${fq_branch%%.fq*}_trimmed.fq.gz
STAR --runThreadN 12 --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn $input --outFileNamePrefix "${fileID}_" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 8000000000 2>&1 | tee -a  ../${fileID}_logs_${dow}.log

echo "cleaning..."
outbam="${fileID}*.sortedByCoord.out.bam"
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
mv $input ../2_read_trimming/

echo "Alignment complete"



