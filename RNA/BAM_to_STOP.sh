#!/bin/bash
set -u

# Script to extract genome coverage across features of interest
# optimised to caluclate 5'-P end frequency adjacent to START or STOP codon from PARE or GMUCT
# only SE

### CONDA environment is installed
# conda create --name ngs_plots
# conda install -n ngs_plots -c bioconda bedtools
# conda install -n ngs_plots r-fields
# conda install -n ngs_plots -c r r-tidyverse
# conda activate ngs_plots

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_STOP.sh <.BAM> <strandedness> <bedfile annotation> <feature name> <distance>"
echo "strandedness = unstranded or  forward"
echo "EXAMPLE: BAM_to_STOP.sh col0_rep1.sorted.bam unstranded Arabidopsis_thaliana.TAIR10.54_stop.bed stop 50"
echo "annotation should be start or stop codons, manually calculated from TAIR10 Ensembl UTR annotations"
exit 1
fi

smp=$1
str=$2
bedfile=$3
out=$4
dis=$5

echo ""
echo "sample = $1"
echo "strand = $2"
echo "bedfile = $3"
echo "feature = $4"
echo "distance = $5"
echo ""

if [[ "$str"  == "unstranded" ]] ; then 

	echo "BAM to bed..."
	bedtools genomecov -bg -5 -ibam $smp > ${smp%%.bam}.5p.bed
	closestBed -D "b" -a ${smp%%.bam}.5p.bed -b $bedfile > ${smp%%.bam}_${out}.5p.bed
	awk -F$'\t' -v a=$dis '$NF<10 && $NF>-a' ${smp%%.bam}_${out}.5p.bed > ${smp%%.bam}_${out}_${dis}bp.5p.bed 

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_stop.r ${smp%%.bam}_${out}_${dis}bp.5p.bed

	echo 'cleaning'
	rm -v ${smp%%.bam}.5p.bed ${smp%%.bam}_${out}.5p.bed

fi

if [[ "$str"  == "forward" ]] ; then

	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand
	
	echo 'stranded BAMs'
	# map to reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# map to foward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	
	echo "BAM to bedgraphs at exons ..."
	# minus strand
	bedtools genomecov -bg -5 -scale -1 -ibam ${smp%%bam}reverse.bam > ${smp%%.bam}.minus.5p.bed
	closestBed -D "b" -a ${smp%%.bam}.minus.5p.bed -b $bedfile > ${smp%%.bam}_${out}.minus.5p.bed
        awk -F$'\t' -v a=$dis '$NF<10 && $NF>-a' ${smp%%.bam}_${out}.minus.5p.bed > ${smp%%.bam}_${out}_${dis}bp.minus.5p.bed

	# plus strand
	bedtools genomecov -bg -5 -ibam ${smp%%bam}forward.bam > ${smp%%.bam}.plus.5p.bed
	closestBed -D "b" -a ${smp%%.bam}.plus.5p.bed -b $bedfile > ${smp%%.bam}_${out}.plus.5p.bed
        awk -F$'\t' -v a=$dis '$NF<10 && $NF>-a' ${smp%%.bam}_${out}.plus.5p.bed > ${smp%%.bam}_${out}_${dis}bp.plus.5p.bed

	echo 'do maths'
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_stop.r ${smp%%.bam}_${out}_${dis}bp.minus.5p.bed
	Rscript /home/dganguly/scripts/RNA/rel_expression_plots_stop.r ${smp%%.bam}_${out}_${dis}bp.plus.5p.bed
	
	echo "Cleaning"
	rm -v  ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%.bam}.minus.5p.bed ${smp%%.bam}_${out}.minus.5p.bed ${smp%%.bam}.plus.5p.bed ${smp%%.bam}_${out}.plus.5p.bed

fi


if [[ "$str"  == "reverse" ]] ; then
	
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand

	echo 'stranded BAMs'
	# reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	
        echo "BAM to bedgraphs at exons ..."
        # minus strand
        bedtools genomecov -bg -5 -scale -1 -ibam ${smp%%bam}forward.bam > ${smp%%.bam}.minus.5p.bed
        closestBed -D "b" -a ${smp%%.bam}.minus.5p.bed -b $bedfile > ${smp%%.bam}_${out}.minus.5p.bed
        awk -F$'\t' -v a=$dis '$NF<10 && $NF>-a' ${smp%%.bam}_${out}.minus.5p.bed > ${smp%%.bam}_${out}_${dis}bp.minus.5p.bed

        # plus strand
        bedtools genomecov -bg -5 -ibam ${smp%%bam}reverse.bam > ${smp%%.bam}.plus.5p.bed
        closestBed -D "b" -a ${smp%%.bam}.plus.5p.bed -b $bedfile > ${smp%%.bam}_${out}.plus.5p.bed
        awk -F$'\t' -v a=$dis '$NF<10 && $NF>-a' ${smp%%.bam}_${out}.plus.5p.bed > ${smp%%.bam}_${out}_${dis}bp.plus.5p.bed

        echo 'do maths'
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots_stop.r ${smp%%.bam}_${out}_${dis}bp.minus.5p.bed
        Rscript /home/dganguly/scripts/RNA/rel_expression_plots_stop.r ${smp%%.bam}_${out}_${dis}bp.plus.5p.bed

        echo "Cleaning"
        rm -v  ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%.bam}.minus.5p.bed ${smp%%.bam}_${out}.minus.5p.bed ${smp%%.bam}.plus.5p.bed ${smp%%.bam}_${out}.plus.5p.bed
	
fi



