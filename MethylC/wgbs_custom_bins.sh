#!/bin/bash
set -eu

# Generate mean methylation levels into custom bins with desired read depth/coverage from per-site BED files

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: wgbs_custom_bins.sh <sample> <file> <genome fasta> <coverage> <bin size>"
echo "EXAMPLE: wgbs_custom_bins.sh col0-r1 bed /home/diep/TAIR10/TAIR10_Chr.all.fasta 15 100"
exit 1
fi

bed=$1
file=$2
fas=$3
cov=$4
bin=$5
window=$(expr $bin - 1)

echo "Weighted methylation in $bed across $bin bp windows with depth >= $cov ..."

cg="${bed}_CG*.${file}"
chg="${bed}_CHG*.${file}"
chh="${bed}_CHH*.${file}"

# use samtools to generate fasta index
samtools faidx $fas

# use awk on index to make genome file
# https://www.biostars.org/p/70795/
awk -v OFS='\t' {'print $1,$2'} ${fas}.fai > temp.genome

# use genome file to make 100bp windows across genome
bedtools makewindows -g temp.genome -w ${window} -s ${bin} | sortBed | awk -F$'\t' ' $1 != "ChrC" && $1 != "ChrM" ' > temp.genome.${bin}bp.sorted.bed

if [ "$file" == "cov" ]; then
# use bedtool intersect and groupBy to get mean methylation levels per bin based on per-site methylation
	echo "Bedtools $cg ..."
	sort -k1,1 -k2,2n $cg | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ),$5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CG_${bin}bp_${cov}cov.bed

	echo "Bedtools $chg ..."
	sort -k1,1 -k2,2n $chg | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ), $5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHG_${bin}bp_${cov}cov.bed

	echo "Bedtools $chh ..."
	sort -k1,1 -k2,2n $chh | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ), $5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHH_${bin}bp_${cov}cov.bed

fi

if [ "$file" == "bed" ]; then
# use bedtool intersect and groupBy to get mean methylation levels per bin based on per-site methylation
        echo "Bedtools $cg ..."
	sort -k1,1 -k2,2n $cg | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 8,9 -o sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($4 / ($4+$5)*100 ),$5 = ($4 + $5)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CG_${bin}bp_${cov}cov.bed

	echo "Bedtools $chg ..."
	sort -k1,1 -k2,2n $chg | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 8,9 -o sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($4 / ($4+$5)*100 ), $5 = ($4 + $5)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHG_${bin}bp_${cov}cov.bed

	echo "Bedtools $chh ..."
	sort -k1,1 -k2,2n $chh | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 8,9 -o sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($4 / ($4+$5)*100 ), $5 = ($4 + $5)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHH_${bin}bp_${cov}cov.bed

fi
	
echo 'cleaning ...'
# CLEAN
rm temp.genome*

