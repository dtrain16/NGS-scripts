#!/bin/bash

# use on sorted BAM file output from bismark alignment to call SNPs using default BS-SNPer settings

## require: 
# perl
# BS-SNPer http://bioinformatics.oxfordjournals.org/content/31/24/4006.long

if [ "$#" -ne 3 ]; then
echo "USAGE: <file> <fasta> <outname>"
echo "EXAMPLE: BS-SNPer.sh alx8-r1.sorted.bam $HOME/TAIR10/TAIR10_Chr.all.fasta alx8-r1"
exit 1
fi

file=$1
fa=$2
out=$3

perl ~/bin/BS-Snper-master/BS-Snper.pl --fa $fa --input $file --output temp.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.15 --minhomfreq 0.85 --minquali 30 --mincover 15 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 > ${out}.SNP.bed 2>${out}_ERR.log

rm temp.out meth.cg meth.chg meth.chh
