#!/bin/bash
set -eu

## Detect peaks from ChIP/RIP-seq relative to a control (input control or RNA-seq) and test for differential enrichment
## Merge replicate all control samples (input, RNA-seq) if replicates are not paired.

### CONDA environment is installed
# conda create --name <name>
# conda install -n <name> -c bioconda bedtools
# conda install -n <name> -c bioconda macs2

if [ "$#" -lt 5 ]; then
echo "Missing required arguments!"
echo "USAGE: MACS_peaks.sh <read layout> <stranded> <input control> <IP> <outname>"
echo "EXAMPLE: MACS_peaks.sh SE/PE un/fr/rf sample1-input_merged.bam sample1-m6_rep1.bam sample1-m6_rep1"
exit 1
fi

#gather input variables
lay=$1
str=$2
input=$3
ip=$4
out=$5

echo $1 $2 $3 $4 $out

if [[ "$lay" == "SE" ]] && [[ "$str"  == "un" ]] ; then

echo "MACS2 Peak calling"

macs2 callpeak --nomodel --extsize 50 -c $input -t $ip -f BAM -g 32542107 -n ${out}.MACS -q 1e-2

fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "fr" ]] ; then

echo "minus bam"

# minus strand
samtools view -@ 4 -f 16 -b $input > ${input%%bam}minus.bam
samtools view -@ 4 -f 16 -b $ip > ${ip%%bam}minus.bam

echo "plus bam"

# plus strand
samtools view -@ 4 -F 16 -b $input > ${input%%bam}plus.bam
samtools view -@ 4 -F 16 -b $ip > ${ip%%bam}plus.bam

echo "Peak calling"

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}plus.bam -t ${ip%%bam}plus.bam -f BAM -g 32542107 -n ${out}_plus.MACS -q 1e-2

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}minus.bam -t ${ip%%bam}minus.bam -f BAM -g 32542107 -n ${out}_minus.MACS -q 1e-2

echo "cleanup"

rm ${input%%bam}plus.bam ${input%%bam}minus.bam ${ip%%bam}plus.bam ${ip%%bam}minus.bam -v

echo "making bedfile of peaks" 

## generate bedfile with strand information using .narrowpeak
awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_minus_peaks.narrowPeak > ${out}_minus_peaks.stranded.bed

awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_plus_peaks.narrowPeak > ${out}_plus_peaks.stranded.bed

## merge and sort by coordinate
cat ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed | sort -k1,1 -k2,2n > ${out}_merged_peaks.strand.bed &

fi


if [[ "$lay" == "SE" ]] && [[ "$str"  == "rf" ]] ; then

echo "minus bam"

# minus strand
samtools view -@ 4 -F 16 -b $input > ${input%%bam}minus.bam
samtools view -@ 4 -F 16 -b $ip > ${ip%%bam}minus.bam

echo "plus bam"

# plus strand
samtools view -@ 4 -f 16 -b $input > ${input%%bam}plus.bam
samtools view -@ 4 -f 16 -b $ip > ${ip%%bam}plus.bam

echo "Peak calling"

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}plus.bam -t ${ip%%bam}plus.bam -f BAM -g 32542107 -n ${out}_plus -q 1e-2

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}minus.bam -t ${ip%%bam}minus.bam -f BAM -g 32542107 -n ${out}_minus -q 1e-2

echo "cleanup"

rm ${input%%bam}plus.bam ${input%%bam}minus.bam ${ip%%bam}plus.bam ${ip%%bam}minus.bam -v

echo "making bedfile of peaks"

## generate bedfile with strand information using .narrowpeak
awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_minus_peaks.narrowPeak > ${out}_minus_peaks.stranded.bed  

awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_plus_peaks.narrowPeak > ${out}_plus_peaks.stranded.bed

## merge and sort by coordinate
cat ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed | sort -k1,1 -k2,2n > ${out}_merged_peaks.strand.bed &

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "un" ]] ; then

echo 'Coming soon!'

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "fr" ]] ; then

echo 'Coming soon!'

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "rf" ]] ; then

echo 'Coming soon!'

fi




