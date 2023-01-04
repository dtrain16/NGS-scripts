#!/bin/bash
set -eu

# Produce files with windowed coverage of RNAseq data (BAM aligned reads) across annotations of interest
# Run in directory with sam converted, sorted, indexed  bam file
# Provide path of genome .fa file to produce windowed genome

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_wigs.sh <.BAM> <genome fasta> <bedfile annotation> <feature name> <window size>"
echo "EXAMPLE: BAM_to_wigs.sh col0_rep1.sorted.bam TAIR10_Chr.all.fasta Araport11_TE.bed TE 100"
exit 1
fi

bam=$1
fas=$2
bedfile=$3
out=$4
size=$5
size_2=$(($size - 1))

echo "Make $size bp genome bed ..."

# use samtools to generate fasta index
samtools faidx $fas
# use awk on index to make genome file
# https://www.biostars.org/p/70795/
awk -v OFS='\t' {'print $1,$2'} ${fas}.fai > temp.genome 
# use genome file to make 100bp windows across genome
bedtools makewindows -g temp.genome -w $size_2 -s $size > temp.genome.${size}bp.bed
sortBed -i temp.genome.${size}bp.bed > temp.genome.${size}bp.sorted.bed

# use bedtools coverage to get coverage across windows from BAM
# MAKE SURE TO USE -sorted FLAG
bedtools coverage -sorted -a temp.genome.${size}bp.bed -b $bam > ${bam%%.sorted*}_${size}bp.bed 

echo 'cleaning ...'
# CLEAN
rm temp.genome*

# sort Bed
sortBed -i ${bam%%.sorted*}_${size}bp.bed > ${bam%%.sorted*}_${size}bp.sorted.bed

echo 'bedtools ...'
# bedtools to desired annotation
closestBed -D "b" -a ${bam%%.sorted*}_${size}bp.sorted.bed -b $bedfile > ${bam%%.sorted*}_${out}_${size}.bed

echo 'subset to +1k/-1k ...'
# awk to subset
awk -F$'\t' '$NF<1000 && $NF>-1000' ${bam%%.sorted*}_${out}_${size}.bed > ${bam%%.sorted*}_${out}_${size}.1k.bed

echo 'final clean ...'
rm ${bam%%.sorted*}_${size}bp.bed ${bam%%.sorted*}_${size}bp.sorted.bed ${bam%%.sorted*}_${out}_${size}.bed

