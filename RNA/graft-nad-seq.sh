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
# conda create --name graft_nad
# conda install -n graft_nad -c bioconda fastqc
# conda install -n graft_nad -c bioconda cutadapt
# conda install -n graft_nad -c bioconda star
# conda install -n graft_nad -c bioconda seqkit 
# conda install -n graft_nad -c conda-forge r-tidyverse
# conda install -n graft_nad conda-forge::parallel ## NOT IMPLEMENTED YET

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

if [[ $fq == *"fq.gz" ]]; then mv ${fq%%.fq*}_fastqc* 1_fastqc; else mv ${fq%%.fastq*}_fastqc* 1_fastqc; fi

echo "Extract branch sequence and trim adapters"
mkdir 1_read_trimming
mv $fq 1_read_trimming
cd 1_read_trimming

## extract reads beginning with branch sequence (GCTTGTTGTG) with flexibility at first and last base
if [[ $fq == *"fq.gz" ]]; 
	then seqkit grep -j 12 -s -r -p "^.CTTGTTGT" $fq -o ${fq%%.fq*}_branch.fq.gz ;
	else seqkit grep -j 12 -s -r -p "^.CTTGTTGT" $fq -o ${fq%%.fastq*}_branch.fq.gz; 
fi

if [[ $fq == *"fastq.gz" ]]; then fq_branch=${fq%%.fastq*}_branch.fq.gz; else fq_branch=${fq%%.fq*}_branch.fq.gz; fi

echo "Trim universal PCR primer sequence from 3' end of read"
## remove universal PCR primer at the 3' end of reads
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" -e 0.2 -m 25 -o "${fq_branch%%.fq*}_3p_trimmed.fq.gz" ${fq_branch} 2>&1 | tee -a ../${fileID}_logs_${dow}.log

## old flags to trim 5' adapter sequences -- obselete
#cutadapt -g "^NCTTGTTGTB" -g "^NCTTGTTGTBB" -g "^NCTTGTTGTBBB" -g  "^NCTTGTTGTBBBG" -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

## R script - trim reads first A at the 5' end of the read, retain read only if A within first 15 bp (10bp branch + flexibility for RT jumping)

echo "Filter and trim reads based on A at 5' end"

## split fq file
val=$(zcat ${fq_branch%%.fq*}_3p_trimmed.fq.gz | wc -l)
val_split=$(expr $bl / 12)
truncated=$(echo "$val_split / 1000000 * 1000000" | bc)
zcat ${fq_branch%%.fq*}_3p_trimmed.fq.gz | split -l $truncated -d --additional-suffix=.fq

## run triming R script on split files in parallel - first numeric argument determines length in which A needs to occur (branch sequence = 10 nts)
parallel -j 12 Rscript ~/scripts/RNA/trim_5p_graft_nad.r 15 {} {}.fq ::: x*

## concatenate output files
#cat *filtered.fq > ${fileID}_filtered_output.fq
cat *trimmed.fq > ${fileID}_trimmed_output.fq
gzip ${fileID}_trimmed_output.fq

# clean up intermediates
rm x*fq
rm *trimmed.fq

## qc filtered and trimmed reads
fastqc -t 12 ${fileID}_trimmed_output.fq.gz 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../
mkdir 0_fastq
mv 1_read_trimming/$fq  0_fastq/
mv 1_read_trimming/$fq_branch 0_fastq/

# read alignment
echo "Align filtered and trimmed reads"

#mkdir 2_align
#mv 1_read_trimming/${fileID}_trimmed_output.fq -t 2_align/
#cd 2_align
#STAR --runThreadN 12 --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 --genomeDir $index --readFilesCommand gunzip -c --readFilesIn ${fileID}_trimmed_output.fq --outFileNamePrefix "${fileID}_" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 8000000000 2>&1 | tee -a  ../${fileID}_logs_${dow}.log
#mv ${fileID}_trimmed_output.fq ../1_read_trimming/

echo "Alignment complete"



