#!/bin/bash
set -u

# Script to summarise read depth from BAM files across features of interest in bedGraph format

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 6 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_bedgraph.sh <.BAM> <layout> <strandedness> <bedfile annotation> <feature name> <distance>"
echo "layout = SE, PE"
echo "strandedness = unstranded, forward, or reverse"
echo "EXAMPLE: BAM_to_bedgraph.sh col0_rep1.sorted.bam SE unstranded Arabidopsis_thaliana.TAIR10.54_exon.bed exon 100"
exit 1
fi

smp=$1
lay=$2
str=$3
bedfile=$4
out=$5
dis=$6

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "bedfile = $4"
echo "feature = $5"
echo "distance = $6"
echo ""

echo "calculate scaling factor"
if [[ "$lay" == "SE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)"); fi
if [[ "$lay" == "PE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -F 260 -c $smp)/2"); fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]] ; then 

	echo "BAM to bed..."
	bedtools genomecov -bg -split -scale $scl -ibam $smp > ${smp%%bam}bed

	echo 'bedtools for coverage across exons...'
	closestBed -D "b" -a ${smp%%bam}bed -b $bedfile > ${smp%%.bam}_${out}.bed

	echo 'subset ...'
	awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.bed > ${smp%%.bam}_${out}.${dis}bp.bed

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.bed

	echo 'cleaning'
	rm -v ${smp%%bam}bed ${smp%%.bam}_${out}.bed ${smp%%.bam}_${out}.${dis}bp.bed

fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "forward" ]] ; then

	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand
	
	echo 'stranded BAMs'
	# map to reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# map to foward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	
	echo "BAM to bedgraphs at exons ..."
	# minus strand
	bedtools genomecov -bg -split -scale -$scl -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bed
	closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bam}_${out}.minus.bed
	awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.minus.bed > ${smp%%.bam}_${out}.${dis}bp.minus.bed

	# plus strand
        bedtools genomecov -bg -split -scale $scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bed
        closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bam}_${out}.plus.bed
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.plus.bed > ${smp%%.bam}_${out}.${dis}bp.plus.bed

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.minus.bed
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.plus.bed
	
	echo "Cleaning"
	rm -v ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.bed ${smp%%.bam}_${out}.minus.bed ${smp%%bam}plus.bed ${smp%%.bam}_${out}.plus.bed

fi


if [[ "$lay" == "SE" ]] && [[ "$str"  == "reverse" ]] ; then
	
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand

	echo 'stranded BAMs'
	# reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	
	echo "BAM to bedgraphs at exons ..."
	# plus strand
	bedtools genomecov -bg -split -scale $scl -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bed
	closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bam}_${out}.plus.bed
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.plus.bed > ${smp%%.bam}_${out}.${dis}bp.plus.bed
        
	# minus strand
        bedtools genomecov -bg -split -scale -$scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bed
        closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bam}_${out}.minus.bed
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.minus.bed > ${smp%%.bam}_${out}.${dis}bp.minus.bed
	
	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.plus.bed
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.minus.bed

	echo "Cleaning ..."	
	rm -v ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.bed ${smp%%.bam}_${out}.minus.bed ${smp%%bam}plus.bed ${smp%%.bam}_${out}.plus.bed

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then
	
	echo "BAM to bed..."
        bedtools genomecov -bg -split -scale $scl -ibam $smp > ${smp%%bam}bed

        echo 'bedtools for coverage across exons...'
        closestBed -D "b" -a ${smp%%bam}bed -b $bedfile > ${smp%%.bam}_${out}.bed

        echo 'subset to +100bp / -100 bp ...'
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.bed > ${smp%%.bam}_${out}.${dis}bp.bed

        Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.bed

	rm ${smp%%bam}bed ${smp%%.bam}_${out}.bed

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "forward" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399
	
	# R1 forward
	samtools view -@ 2 -f 99 -b $smp > ${smp%%bam}R1F.bam
	# R2 reverse
	samtools view -@ 2 -f 147 -b $smp > ${smp%%bam}R2R.bam
	# FORWARD R1 read pairs
	samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

	# R1 reverse
	samtools view -@ 2 -f 83 -b $smp > ${smp%%bam}R1R.bam
	# R2 forward
	samtools view -@ 2 -f 163 -b $smp > ${smp%%bam}R2F.bam
	# REVERSE R1 read pairs
	samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	rm ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	echo "BAM to bedgraph at feature..."
	# minus strand
	bedtools genomecov -bg -split -scale -$scl -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bed
	closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bam}_${out}.minus.bed	
	awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.minus.bed > ${smp%%.bam}_${out}.${dis}bp.minus.bed

	# plus strand
	bedtools genomecov -bg -split -scale $scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bed
	closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bam}_${out}.plus.bed
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.plus.bed > ${smp%%.bam}_${out}.${dis}bp.plus.bed
	
	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.minus.bed
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.plus.bed

	rm ${smp%%bam}forward.bam ${smp%%bam}reverse.bam ${smp%%bam}minus.bed ${smp%%.bam}_${out}.minus.bed ${smp%%bam}plus.bed  ${smp%%.bam}_${out}.plus.bed

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "reverse" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399

	# R1 forward
	samtools view -@ 2 -f 99 -b $smp > ${smp%%bam}R1F.bam
	# R2 reverse
	samtools view -@ 2 -f 147 -b $smp > ${smp%%bam}R2R.bam
	# FORWARD R1 read pairs
	samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

	# R1 reverse
	samtools view -@ 2 -f 83 -b $smp > ${smp%%bam}R1R.bam
	# R2 forward
	samtools view -@ 2 -f 163 -b $smp > ${smp%%bam}R2F.bam
	# REVERSE R1 read pairs
	samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	rm ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	echo "BAM to stranded bedgraph ..."
	# plus strand
	bedtools genomecov -bg -split -scale $scl -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bed
	closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bam}_${out}.plus.bed
	awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.plus.bed > ${smp%%.bam}_${out}.${dis}bp.plus.bed	

	# minus strand
	bedtools genomecov -bg -split -scale -$scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bed
	closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bam}_${out}.minus.bed
        awk -F$'\t' -v a=$dis '$NF<a && $NF>-a' ${smp%%.bam}_${out}.minus.bed > ${smp%%.bam}_${out}.${dis}bp.minus.bed

	echo 'do maths'
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.minus.bed
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bam}_${out}.${dis}bp.plus.bed
	
	rm ${smp%%.bam}_${out}.minus.bed ${smp%%bam}minus.bed ${smp%%.bam}_${out}.plus.bed ${smp%%bam}plus.bed ${smp%%bam}reverse.bam ${smp%%bam}forward.bam

fi


