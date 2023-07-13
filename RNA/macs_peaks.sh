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

echo "making bedfile of peaks" 

## generate bedfile with strand information using .narrowpeak
awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_minus_peaks.narrowPeak > ${out}_minus_peaks.stranded.bed

awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_plus_peaks.narrowPeak > ${out}_plus_peaks.stranded.bed

## merge and sort by coordinate
cat ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed | sort -k1,1 -k2,2n > ${out}_merged_peaks.strand.bed &

echo "cleanup"

rm -v ${input%%bam}plus.bam ${input%%bam}minus.bam ${ip%%bam}plus.bam ${ip%%bam}minus.bam
rm -v ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed

echo "Done!"

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

echo "making bedfile of peaks"

## generate bedfile with strand information using .narrowpeak
awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_minus_peaks.narrowPeak > ${out}_minus_peaks.stranded.bed  

awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_plus_peaks.narrowPeak > ${out}_plus_peaks.stranded.bed

## merge and sort by coordinate
cat ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed | sort -k1,1 -k2,2n > ${out}_merged_peaks.strand.bed &

echo "cleanup"

rm -v ${input%%bam}plus.bam ${input%%bam}minus.bam ${ip%%bam}plus.bam ${ip%%bam}minus.bam
rm -v ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed

echo "Done!"

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "un" ]] ; then

echo 'Coming soon!'

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "fr" ]] ; then

echo 'Coming soon!'

# R1 forward
samtools view -@ 2 -f 99 -b $smp > ${smp%%bam}R1F.bam
# R2 reverse
samtools view -@ 2 -f 147 -b $smp > ${smp%%bam}R2R.bam
# FORWARD R1 read pairs
samtools merge -f ${smp%%bam}plus.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

# R1 reverse
samtools view -@ 2 -f 83 -b $smp > ${smp%%bam}R1R.bam
# R2 forward
samtools view -@ 2 -f 163 -b $smp > ${smp%%bam}R2F.bam
# REVERSE R1 read pairs
samtools merge -f ${smp%%bam}minus.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

rm ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "rf" ]] ; then

echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
# http://seqanswers.com/forums/showthread.php?t=29399

echo "input samples"

# R1 forward
samtools view -@ 2 -f 99 -b $input > ${input%%bam}R1F.bam
# R2 reverse
samtools view -@ 2 -f 147 -b $input > ${input%%bam}R2R.bam
# FORWARD R1 read pairs
samtools merge -f ${input%%bam}minus.bam ${input%%bam}R1F.bam ${input%%bam}R2R.bam

# R1 reverse
samtools view -@ 2 -f 83 -b $input > ${input%%bam}R1R.bam
# R2 forward
samtools view -@ 2 -f 163 -b $input > ${input%%bam}R2F.bam
# REVERSE R1 read pairs
samtools merge -f ${input%%bam}plus.bam ${input%%bam}R1R.bam ${input%%bam}R2F.bam

rm ${input%%bam}R1F.bam ${input%%bam}R2R.bam ${input%%bam}R1R.bam ${input%%bam}R2F.bam

echo "IP samples"

# R1 forward
samtools view -@ 2 -f 99 -b $ip > ${ip%%bam}R1F.bam
# R2 reverse
samtools view -@ 2 -f 147 -b $ip > ${ip%%bam}R2R.bam
# FORWARD R1 read pairs
samtools merge -f ${ip%%bam}minus.bam ${ip%%bam}R1F.bam ${ip%%bam}R2R.bam

# R1 reverse
samtools view -@ 2 -f 83 -b $ip > ${ip%%bam}R1R.bam
# R2 forward
samtools view -@ 2 -f 163 -b $ip > ${ip%%bam}R2F.bam
# REVERSE R1 read pairs
samtools merge -f ${ip%%bam}plus.bam ${ip%%bam}R1R.bam ${ip%%bam}R2F.bam

rm ${ip%%bam}R1F.bam ${ip%%bam}R2R.bam ${ip%%bam}R1R.bam ${ip%%bam}R2F.bam

echo 'peak calling'

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}plus.bam -t ${ip%%bam}plus.bam -f BAM -g 32542107 -n ${out}_plus -q 1e-2

macs2 callpeak --nomodel --extsize 50 -c ${input%%bam}minus.bam -t ${ip%%bam}minus.bam -f BAM -g 32542107 -n ${out}_minus -q 1e-2

## generate bedfile with strand information using .narrowpeak
awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "-"}' ${out}_minus_peaks.narrowPeak > ${out}_minus_peaks.stranded.bed

awk 'BEGIN {OFS="\t";FS="\t"} $ awk {print $1,$2,$3,$4,$5, "+"}' ${out}_plus_peaks.narrowPeak > ${out}_plus_peaks.stranded.bed

## merge and sort by coordinate
cat ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed | sort -k1,1 -k2,2n > ${out}_merged_peaks.strand.bed &

echo "cleanup"

rm -v ${input%%bam}plus.bam ${input%%bam}minus.bam ${ip%%bam}plus.bam ${ip%%bam}minus.bam ${out}_minus_peaks.stranded.bed ${out}_plus_peaks.stranded.bed

echo "Done!"


fi


