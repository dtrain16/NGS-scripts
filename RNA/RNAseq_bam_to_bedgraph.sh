#!/bin/bash
set -e 
set -u

# Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format
# Then produce bigWigs files for viewing delight
# Run in directory with sam converted, sorted, indexed  bam file
# Ensure subread indexed genome & chromosome sizes are prepared eg TAIR10/TAIR10_Chr.all.fasta.len
# samtools faidx TAIR10_Chr.all.fasta | cut -f1,2 TAIR10_Chr.all.fasta.fai > TAIR10_Chr.all.fasta.len

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_bam_to_bedgraph.sh <sample name> <SE, PE> <unstranded, stranded, rev_stranded> <chr_sizes>"
echo "EXAMPLE: RNAseq_bam_to_bedgraph.sh col0-r1.bam PE unstranded TAIR10_Chr.all.fasta.len"
exit 1
fi

smp=$1
lay=$2
str=$3
chrc_sizes=$4

# file for length of all 7 chromosomes
# chrc_sizes=$HOME/TAIR10/TAIR10_Chr.all.fasta.len

echo ""
echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "chr_size = $4"
echo ""
echo ""
echo "Produce $lay $str bigWig file(s) from $smp ..."
echo ""

if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]] ; then 

	echo "BAM to bedgraph ..."
	# unstranded bedgraph
	bedtools genomecov -bga -split -ibam $smp -g $chrc_sizes > ${smp%%bam}bg

	# bg to bigWig
	echo "bigWig ..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}bg ${chrc_sizes} ${smp%%bam}bigWig

fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "stranded" ]] ; then
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand
	
	# reverse strand
	samtools view -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -F 16 -b $smp > ${smp%%bam}forward.bam
	echo "BAM to stranded bedgraphs ..."
	# reverse/minus bg
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}reverse.bam -g $chrc_sizes > ${smp%%bam}minus.bg
	# forward/plus bg
	bedtools genomecov -bga -split -ibam ${smp%%bam}forward.bam -g $chrc_sizes > ${smp%%bam}plus.bg
	
	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes} ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig
	
	rm ${smp%%.bam}*forward*bam -v
	rm ${smp%%.bam}*reverse*bam -v

fi


if [[ "$lay" == "SE" ]] && [[ "$str"  == "rev_stranded" ]] ; then
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand

	# reverse strand
	samtools view -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -F 16 -b $smp > ${smp%%bam}forward.bam
	echo "BAM to stranded bedgraphs ..."
	# reverse/plus bg
	bedtools genomecov -bga -split -ibam ${smp%%bam}reverse.bam -g $chrc_sizes > ${smp%%bam}plus.bg
	# forward/minus bg
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}forward.bam -g $chrc_sizes > ${smp%%bam}minus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig
	rm ${smp%%.bam}*forward*bam -v
	rm ${smp%%.bam}*reverse*bam -v

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then
	# need sorted bam
	echo "sort by position"
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"
	
	echo "BAM to bedgraph ..."
	# unstranded bedgraph
	bedtools genomecov -bga -split -ibam $smp > ${smp%%bam}bg

	# bg to bigWig
	echo "bigWig ..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}bg ${chrc_sizes} ${smp%%bam}bigWig

	rm $smp ${smp%%bam}bg -v

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "stranded" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399
	
	# need sorted bam
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"

	# R1 forward
	samtools view -f 99 -b $smp > ${smp%%bam}R1F.bam
	# R2 reverse
	samtools view -f 147 -b $smp > ${smp%%bam}R2R.bam
	# FORWARD R1 read pairs
	samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

	# R1 reverse
	samtools view -f 83 -b $smp > ${smp%%bam}R1R.bam
	# R2 forward
	samtools view -f 163 -b $smp > ${smp%%bam}R2F.bam
	# REVERSE R1 read pairs
	samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	echo "BAM to stranded bedgraph ..."
	# minus strand
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bg
	# plus strand
	bedtools genomecov -bga -split -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig

	rm ${smp%%.bam}*R1*bam ${smp%%.bam}*R2*bam -v
	rm ${smp%%.bam}*forward*bam ${smp%%.bam}*reverse*bam $smp -v

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "rev_stranded" ]] ; then

	echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
	# http://seqanswers.com/forums/showthread.php?t=29399

	# need sorted bam
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"

	# R1 forward
	samtools view -f 99 -b $smp > ${smp%%bam}R1F.bam
	# R2 reverse
	samtools view -f 147 -b $smp > ${smp%%bam}R2R.bam
	# FORWARD R1 read pairs
	samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

	# R1 reverse
	samtools view -f 83 -b $smp > ${smp%%bam}R1R.bam
	# R2 forward
	samtools view -f 163 -b $smp > ${smp%%bam}R2F.bam
	# REVERSE R1 read pairs
	samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

	echo "BAM to stranded bedgraph ..."
	# plus strand
	bedtools genomecov -bga -split -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bg
	# minus strand
	bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig

	rm ${smp%%.bam}*R1*bam ${smp%%.bam}*R2*bam -v
	rm ${smp%%.bam}*forward*bam ${smp%%.bam}*reverse*bam -v
	rm $smp -v

fi

# clean up tmps
rm ${smp%%.sorted.bam}*bg -v
