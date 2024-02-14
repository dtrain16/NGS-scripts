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

if [[ "$lay" == "SE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -c $smp)"); fi
if [[ "$lay" == "PE" ]] ; then scl=$(bc <<< "scale=6;1000000/$(samtools view -f 1 -c $smp)"); fi

echo ""
echo "sample = $1"
echo "chr_size = $2"
echo ""
echo "Produce bigWig file(s) for 5p read ends from $smp ..."
echo ""


echo "BAM to bedgraph ..."
# unstranded bedgraph counting only 5p read end
bedtools genomecov -bga -5 -scale $scl -ibam $smp > ${smp%%bam}5p.bg

# bg to bigWig
echo "bigWig ..."
$HOME/bin/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${smp%%bam}5p.bg ${chrc_sizes} ${smp%%bam}5p.bigWig

rm ${smp%%bam}5p.bg

