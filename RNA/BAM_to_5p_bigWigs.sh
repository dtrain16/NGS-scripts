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

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: BAM_to_5p_bigWig.sh <BAM> <layout> <chr_sizes>"
echo "EXAMPLE: BAM_to_5p_bigWig.sh col0-r1.bam SE TAIR10_Chr.all.fasta.len"
exit 1
fi

smp=$1
lay=$2
chrc_sizes=$3

echo ""
echo "sample = $1"
echo "layout = $2"
echo "chr_size = $3"
echo ""
echo "Produce bigWig file(s) for 5p read ends from $smp ..."
echo ""


if [[ "$lay" == "SE" ]] ; then 

reads=$(samtools view -F 260 -c $smp)
scaling_factor=$(bc <<< "scale=6;1000000/$reads")

echo "BAM to bedgraph ..."
# unstranded bedgraph of 5' read end coverage unscaled
bedtools genomecov -bga -5 -ibam $smp > ${smp%%bam}5p.raw.bg

# unstranded bedgraph of 5' read end coverage scaled to RPM
bedtools genomecov -bga -5 -scale $scaling_factor -ibam $smp > ${smp%%bam}5p.rpm.bg

fi

if [[ "$lay" == "PE" ]] ; then

# unstranded bedgraph of 5' read end coverage unscaled
bedtools genomecov -bga -5 -ibam $smp > ${smp%%bam}5p.raw.bg

reads=$(samtools view -F 260 -c $smp)
frags=$(expr $reads / 2)
scaling_factor=$(bc <<< "scale=6;1000000/$frags")

echo "BAM to bedgraph ..."
# unstraned bedgraph of 5' read end coverage scaled to RPM
bedtools genomecov -bga -5 -scale $scaling_factor -ibam $smp > ${smp%%bam}5p.rpm.bg

fi

# convert bedgraph to bigWig
echo "bigWig ..."
$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}5p.raw.bg ${chrc_sizes} ${smp%%bam}5p.raw.bigWig

$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}5p.rpm.bg ${chrc_sizes} ${smp%%bam}5p.rpm.bigWig

rm ${smp%%bam}5p.raw.bg ${smp%%bam}5p.rpm.bg


