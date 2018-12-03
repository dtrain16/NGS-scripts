#!/bin/bash
set -eu
bam=$1

java -jar ~/bin/picard.jar CleanSam I=$bam O=test.bam
java -jar ~/bin/picard.jar ValidateSamFile I=test.bam IGNORE_WARNINGS=true MODE=VERBOSE

# AddOrReplaceReadGroups
# FixMateInformation
## Manually remove reads with errors
# samtools view -h test.bam | grep -v 'D00775:83:CC65GANXX:4:2304:19617:7489' | samtools view -b > test2.bam

samtools sort -@ 6 test.bam -o test.sorted.bam

java -jar ~/bin/picard.jar CollectInsertSizeMetrics I=test.sorted.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf

rm test.bam
rm test.sorted.bam

echo "DONE"
