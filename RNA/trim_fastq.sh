#!/bin/bash

## hard trim FASTQs sequences e.g. for PARE-seq libraries
## https://wikis.utexas.edu/display/CoreNGSTools/Pre-processing+raw+sequences#Preprocessingrawsequences-FASTXToolkit

### CONDA environment is installed
# conda create --name <name>
# conda install -c bioconda fastx_toolkit


if [ "$#" -lt 2 ]; then
echo "Missing required arguments!"
echo "USAGE: trim_fastq.sh <fastq> <length>"
echo "EXAMPLE: trim_fastq.sh sample.fastq.gz 15"
exit 1
fi

fq=$1
NSEQS=$2

if [[ $fq != *.gz ]];then fastx_trimmer -z -l $NSEQS -i $fq -o ${fq%%.fastq}.${NSEQS}bp.fastq.gz; fi

if [[ $fq == *.gz ]];then zcat $fq | fastx_trimmer -z -l $NSEQS -o ${fq%%.fastq.gz}.${NSEQS}bp.fastq.gz; fi


