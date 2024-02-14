#!/bin/bash
set -eu

# Produce scaled coverage data (in RPM) from BAM files in bedgraph format.
# Then produce bigWigs files for viewing delight
# Run in directory with sam converted, sorted, indexed  bam file
# Ensure genome index genome & chromosome sizes are prepared:
# samtools faidx TAIR10_Chr.all.fasta | cut -f1,2 TAIR10_Chr.all.fasta.fai > TAIR10_Chr.all.fasta.len
# Make sure you have kentUtils: https://github.com/ENCODE-DCC/kentUtils
# git clone git://github.com/ENCODE-DCC/kentUtils.git

### CONDA environment is installed
# conda create --name Bedtools
# conda install -n Bedtools -c bioconda bedtools

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_bigWig.sh <BAM> <SE, PE> <unstranded, stranded, rev_stranded> <chr_sizes>"
echo "EXAMPLE: BAM_to_bigWig.sh col0-r1.bam PE unstranded TAIR10_Chr.all.fasta.len"
exit 1
fi

smp=$1
lay=$2
str=$3
chrc_sizes=$4

if [[ "$lay" == "SE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -c $smp)"); fi
if [[ "$lay" == "PE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -f 1 -c $smp)"); fi

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "chr_size = $4"
echo ""
echo "Produce $lay $str bigWig file(s) from $smp ..."
echo ""

if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]] ; then 

	echo "BAM to bedgraph ..."
	# unstranded bedgraph
	bedtools genomecov -bga -split -scale $scl -ibam $smp > ${smp%%bam}bg

	# bg to bigWig
	echo "bigWigigWig ..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}bg ${chrc_sizes} ${smp%%bam}bigWig

	rm ${smp%%bam}bg

fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "stranded" ]] ; then
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand
	
	# reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	echo "BAM to stranded bedgraphs ..."
	# reverse/minus bg
	bedtools genomecov -bga -split -scale -${scl} -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bg
	# forward/plus bg
	bedtools genomecov -bga -split -scale $scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bg
	
	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes} ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig
	
	rm ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.bg ${smp%%bam}plus.bg

fi


if [[ "$lay" == "SE" ]] && [[ "$str"  == "rev_stranded" ]] ; then
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand

	# reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	Samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	echo "BAM to stranded bedgraphs ..."
	# reverse/plus bg
	bedtools genomecov -bga -split -scale $scl -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bg
	# forward/minus bg
	bedtools genomecov -bga -split -scale -${scl} -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig
	
	rm ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}plus.bg ${smp%%bam}minus.bg

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then
	# need sorted bam
	echo "sort by position"
	samtools sort -@ 4 ${smp} -o ${smp%%bam}sorted.bam
	smp="${smp%%bam}sorted.bam"
	
	echo "BAM to bedgraph ..."
	# unstranded bedgraph
	bedtools genomecov -bga -split -scale $scl -ibam $smp > ${smp%%bam}bg

	# bg to bigWig
	echo "bigWig ..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}bg ${chrc_sizes} ${smp%%bam}bigWig

	rm ${smp%%bam}bg $smp

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "stranded" ]] ; then

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
	bedtools genomecov -bga -split -scale -${scl} -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.bg
	# plus strand
	bedtools genomecov -bga -split -scale ${scl} -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig

	rm $smp ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}forward.bam ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam ${smp%%bam}minus.bg ${smp%%bam}plus.bg

fi

if [[ "$lay" == "PE" ]] && [[ "$str"  == "rev_stranded" ]] ; then

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
	bedtools genomecov -bga -split -scale ${scl} -ibam ${smp%%bam}reverse.bam > ${smp%%bam}plus.bg
	# minus strand
	bedtools genomecov -bga -split -scale -${scl} -ibam ${smp%%bam}forward.bam > ${smp%%bam}minus.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig

	rm ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam ${smp%%bam}plus.bg ${smp%%bam}minus.bg

fi

