#!/bin/bash
set -eu

## Detect peaks from m6A-RNA-seq and test for differential enrichment

### CONDA environment is installed
# conda create --name <name>
# conda install -n <name> -c bioconda bedtools
# conda install -n <name> -c bioconda macs2

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: m6a_peaks.sh <SE/PE> <fastq R1> <R2> </path/to/index> <fileID>"
echo "EXAMPLE: m6a_peaks SE sample.fastq /home/dganguly/ref_seqs/STAR/TAIR10/GenomeDir sample_rep1"
exit 1
fi


#gather input variables
type=$1
fq=$2;
index=$3; #path to subread indexed reference genome
fileID=$4;



macs2 callpeak --nomodel --extsize 50 -c /Data05/prallw/redo_mRNAseq/STAR/10B_1_C_4hr_pairedAligned.sortedByCoord.out.plus.bam -t /Data05/prallw/redo_m6Aseq/STAR/10B_m6_1_C_0hr_pairedAligned.sortedByCoord.out.plus.bam -f BAM -g 32542107 -n /Data05/prallw/redo_m6Aseq/MACS/10B_m6_1_C_0hr_MACS.plus -q 1e-2




