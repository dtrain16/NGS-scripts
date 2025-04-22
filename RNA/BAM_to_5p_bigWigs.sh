set -eu

# Produce 5p end coverage data from BAM files from GMUCT or PARE-seq in bedgraph format
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
echo "USAGE: BAM_to_5p_bigWig.sh <BAM> <layout> <strand> <chr_sizes>"
echo "EXAMPLE: BAM_to_5p_bigWig.sh col0-r1.bam SE unstranded,stranded TAIR10_Chr.all.fasta.len"
exit 1
fi

smp=$1
lay=$2
str=$3
chrc_sizes=$4

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "chr_size = $4"
echo ""
echo "Produce bigWig file(s) for 5p read ends from $smp ..."
echo ""


if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]]; then

	reads=$(samtools view -F 260 -c $smp)
	scaling_factor=$(bc <<< "scale=6;1000000/$reads")

	echo "BAM to bedgraph ..."
	# unstranded bedgraph of 5' read end coverage scaled to RPM
	bedtools genomecov -bga -5 -scale $scaling_factor -ibam $smp > ${smp%%bam}5p.bg

	# convert bedgraph to bigWig
	echo "bigWig ..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}5p.bg ${chrc_sizes} ${smp%%bam}5p.bigWig


fi

if [[ "$lay" == "SE" ]] && [[ "$str"  == "stranded" ]] ; then
	# https://www.biostars.org/p/179035/
	# extract reads from + and - strand
	
	reads=$(samtools view -F 260 -c $smp)
	scaling_factor=$(bc <<< "scale=6;1000000/$reads")

	# reverse strand
	samtools view -@ 2 -f 16 -b $smp > ${smp%%bam}reverse.bam
	# forward strand
	samtools view -@ 2 -F 16 -b $smp > ${smp%%bam}forward.bam
	
	echo "BAM to stranded bedgraphs ..."
	# reverse/minus bg
	bedtools genomecov -bga -5 -scale -${scl} -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.5p.bg
	# forward/plus bg
	bedtools genomecov -bga -5 -scale $scl -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.5p.bg
	
	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.5p.bg ${chrc_sizes} ${smp%%bam}plus.5p.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.5p.bg ${chrc_sizes} ${smp%%bam}minus.5p.bigWig
	
	rm ${smp%%bam}reverse.bam ${smp%%bam}forward.bam ${smp%%bam}minus.5p.bg ${smp%%bam}plus.5p.bg

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then

reads=$(samtools view -F 260 -c $smp)
frags=$(expr $reads / 2)
scaling_factor=$(bc <<< "scale=6;1000000/$frags")

echo "BAM to bedgraph ..."
# unstraned bedgraph of 5' read end coverage scaled to RPM
bedtools genomecov -bga -5 -scale $scaling_factor -ibam $smp > ${smp%%bam}5p.bg

# convert bedgraph to bigWig
echo "bigWig ..."
$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}5p.bg ${chrc_sizes} ${smp%%bam}5p.bigWig

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "stranded" ]] ; then

	reads=$(samtools view -F 260 -c $smp)
	frags=$(expr $reads / 2)
	scaling_factor=$(bc <<< "scale=6;1000000/$frags")


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
	
	rm $smp ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam ${smp%%bam}forward.bam ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam	

	echo "BAM to stranded bedgraph ..."
	# minus strand
	bedtools genomecov -bga -5 -scale -${scl} -ibam ${smp%%bam}reverse.bam > ${smp%%bam}minus.5p.bg
	# plus strand
	bedtools genomecov -bga -5 -scale ${scl} -ibam ${smp%%bam}forward.bam > ${smp%%bam}plus.5p.bg

	echo "bigWigs..."
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}plus.5p.bg ${chrc_sizes}  ${smp%%bam}plus.5p.bigWig
	$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}minus.5p.bg ${chrc_sizes} ${smp%%bam}minus.5p.bigWig

	rm ${smp%%bam}minus.5p.bg ${smp%%bam}plus.5p.bg

fi


