#!/bin/bash
set -eu

# Script to summarise read depth from BAM files across features of interest in bedGraph format

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_bedgraph.sh <.BAM> <layout> <strandedness> <bedfile annotation> <feature name>"
echo "layout = SE, PE"
echo "strandedness = unstranded, forward, or reverse"
echo "EXAMPLE: BAM_to_bedgraph.sh col0_rep1.sorted.bam SE unstranded Arabidopsis_thaliana.TAIR10.54_exon.bed exon"
exit 1
fi

smp=$1
lay=$2
str=$3
bedfile=$4
out=$5

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "bedfile = $4"
echo "feature = $5"
echo ""

if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]] ; then 

	echo "BAM to bed..."
	bedtools genomecov -bga -split -ibam $smp > ${smp%%bam}bed

	echo 'bedtools for coverage across exons...'
	closestBed -D "b" -a ${smp%%bam}bed -b $bedfile > ${smp%%.bed*}_${out}.bed

	echo 'subset ...'
	awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.bed > ${smp%%.bed*}_${out}.10bp.bed

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc_v1.r ${smp%%.bed*}_${out}.10bp.bed

	echo 'cleaning'
	rm -v ${smp%%bam}bed ${smp%%.bed*}_${out}.bed

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
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bed
	closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bed*}_${out}.minus.bed
	awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.minus.bed > ${smp%%.bed*}_${out}.10bp.minus.bed

	# plus strand
        bedtools genomecov -bga -split -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bed
        closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bed*}_${out}.plus.bed
        awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.plus.bed > ${smp%%.bed*}_${out}.10bp.plus.bed

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc_v1.r ${smp%%.bed*}_${out}.10bp.minus.bed
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc_v1.r ${smp%%.bed*}_${out}.10bp.plus.bed
	
	echo "Cleaning"
	rm -v ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.bed ${smp%%.bed*}_${out}.minus.bed ${smp%%bam}plus.bed ${smp%%.bed*}_${out}.plus.bed

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
	bedtools genomecov -bga -split -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bed
	closestBed -D "b" -a ${smp%%bam}plus.bed -b $bedfile > ${smp%%.bed*}_${out}.plus.bed
        awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.plus.bed > ${smp%%.bed*}_${out}.10bp.plus.bed
        
	# minus strand
        bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bed
        closestBed -D "b" -a ${smp%%bam}minus.bed -b $bedfile > ${smp%%.bed*}_${out}.minus.bed
        awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.minus.bed > ${smp%%.bed*}_${out}.10bp.minus.bed
	
	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc_v1.r ${smp%%.bed*}_${out}.10bp.plus.bed
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots_ejc_v1.r ${smp%%.bed*}_${out}.10bp.minus.bed

	echo "Cleaning ..."	
	rm -v ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.bed ${smp%%.bed*}_${out}.minus.bed ${smp%%bam}plus.bed ${smp%%.bed*}_${out}.plus.bed

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then
	
	echo "PE not developed properly!!!"
	
	echo "BAM to bed..."
        bedtools genomecov -bga -split -ibam $smp > ${smp%%.bam*}.bed

        echo 'bedtools for coverage across exons...'
        closestBed -D "b" -a ${smp%%.bam*}.bed -b $bedfile > ${smp%%.bed*}_${out}.bed

        echo 'subset to + 50bp / -50 bp ...'
        awk -F$'\t' '$NF<10 && $NF>-10' ${smp%%.bed*}_${out}.bed > ${smp%%.bed*}_${out}.10bp.bed

        Rscript /home/dganguly/scripts/RNA/rel_expression_plots.r ${smp%%.bed*}_${out}.10bp.bed

	

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "forward" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399
	
	# need sorted bam
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"

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

	echo "BAM to stranded bedgraph ..."
	# minus strand
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bg
	# plus strand
	bedtools genomecov -bga -split -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bg

	rm ${smp%%.bam}*R1*bam ${smp%%.bam}*R2*bam -v
	rm ${smp%%.bam}*forward*bam ${smp%%.bam}*reverse*bam $smp -v

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "reverse" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399

	# need sorted bam
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"

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

	echo "BAM to stranded bedgraph ..."
	# plus strand
	bedtools genomecov -bga -split -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bg
	# minus strand
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bg

	rm ${smp%%.bam}*R1*bam ${smp%%.bam}*R2*bam -v
	rm ${smp%%.bam}*forward*bam ${smp%%.bam}*reverse*bam -v
	rm $smp -v

fi


