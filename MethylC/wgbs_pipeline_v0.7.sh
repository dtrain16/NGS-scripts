#!/bin/bash
set -u

###################
# This script it designed to take a single end fastq file and process it through the 
# bismark aligner call methylated cytosines, and develop per-c bed files and 100bp 
# window wig files for CG, CHG, and CHH methylation levels
# Compatible with Bismark v0.19
# USES BOWTIE2 FOR ALIGNMENT
# Genome indexing
# Bowtie2: bismark_genome_preparation --bowtie2 /path/to/genome
###################

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: wgbs_pipelinev0.6.sh <SE/PE> <R1> <R2> <path to bismark genome> <fileID for output>"
exit 1
fi

######################################
# SINGLE END
######################################

#confirm single-end
if [ "$1" == "SE" ];then

#require arguments
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: wgbs_pipelinev0.6.sh SE <R1> <path to bismark genome> <fileID for output>"
exit 1
fi

#gather input variables
type=$1; #identifying paired end or single end mode
fq_file=$2; #the input fastq file
genome_path=$3; #the path to the genome to be used (bismark genome prepped)
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing Bismark single-end alignment with the following parameters:"
echo "Type: $type"
echo "Input File: $fq_file"
echo "Path to bismark genome folder: $genome_path"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo ""
echo "Full logfile of steps: ${fileID}_logs_${dow}.log"
echo "##################"

#develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

#fastqc
mkdir 1_fastqc
fastqc -t 4 $fq_file 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq_file%%.fastq*}_fastqc* 1_fastqc

## trim_galore!
mkdir 2_trimgalore
cd 2_trimgalore
trim_galore ../$fq_file 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../

#fastqc_again
mkdir 3_trimmed_fastqc
fastqc -t 4 2_trimgalore/${fq_file%%.fastq*}_trimmed.fq* 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq_file%%.fastq*}_trimmed_fastqc* 3_trimmed_fastqc

mkdir 0_rawfastq
mv $fq_file 0_rawfastq

## bismark to BAM
mkdir 4_bismark_alignment
cd 4_bismark_alignment

bismark --multicore 2 $genome_path ../2_trimgalore/${fq_file%%.fastq*}_trimmed.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

## deduplicate BAMs (e.g. PCR duplicates)
deduplicate_bismark --bam -s ${fq_file%%.fastq*}_trimmed_bismark_bt2.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# sort BAM
samtools sort ${fq_file%%.fastq*}_trimmed_bismark_bt2.deduplicated.bam \
	-o ${fq_file%%.fastq*}_bt2.deduplicated.sorted.bam

samtools index ${fq_file%%.fastq*}_bt2.deduplicated.sorted.bam

# remove intermediate BAM
rm -v ${fq_file%%.fastq*}_trimmed_bismark_bt2.bam
rm -v ${fq_file%%.fastq*}_trimmed_bismark_bt2.deduplicated.bam

## methylation extraction
bismark_methylation_extractor --comprehensive --report --multicore 3 --buffer_size 8G \
	-s ${fq_file%%.fastq*}_bt2.deduplicated.sorted.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

## bedgraph creation
bismark2bedGraph --CX CpG*txt -o ${fileID}_CG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHG*txt -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHH*txt -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

rm -v C*txt
rm -v *bed.gz

cd ../

mkdir 5_output_files
mv 4_bismark_alignment/*.cov.gz 5_output_files

echo "single C bed file with depth > 3"
## depth filter on single Cs
cd 5_output_files

zcat ${fileID}_CG.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

zcat ${fileID}_CHG.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

zcat ${fileID}_CHH.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../

echo "#####################"
echo "providing pipeline metrics to wgbs pipeline logfile..."
echo "#####################"

#get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of alignments with a unique best hit from the different alignments:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequences with no alignments under any condition:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequences did not map uniquely:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed_bismark_bt2_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

if [[ $fq_file == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4))
fi

if [[ $fq_file != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq_file%%.fastq*}_trimmed.fq*)
	flt_reads=$(($flt_reads / 4))
fi

## add it to the full pipeline logfile
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path/}\t${type}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path/}\t${type}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log

fi

######################################
# PAIRED END
######################################

if [ "$1" == "PE" ];then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: wgbs_pipelinev0.6.sh PE <R1> <R2> <path to bismark genome folder> <fileID for output files>"
exit 1
fi
#gather input variables
type=$1; #identifying paired end or single end mode
fq_file1=$2; #R1 reads
fq_file2=$3; #R2 reads
genome_path=$4; #the path to the genome to be used (bismark genome prepped)
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing Bismark paired-end alignment with the following parameters:"
echo "Type: $type"
echo "Input File: $fq_file1"
echo "Path to bismark genome folder: $genome_path"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo ""
echo "Full logfile of steps: ${fileID}_logs_${dow}.log"
echo "##################"

## develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file1 ${fileID}_wgbspipeline_${dow}
mv $fq_file2 ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

## fastqc
mkdir 1_fastqc
fastqc -t 4 $fq_file1 $fq_file2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq_file1%%.fastq*}_fastqc* 1_fastqc
mv ${fq_file2%%.fastq*}_fastqc* 1_fastqc

## trim_galore!
mkdir 2_trimgalore
cd 2_trimgalore
trim_galore --paired ../$fq_file1 ../$fq_file2 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../

## fastqc_again
mkdir 3_trimmed_fastqc
fastqc -t 4 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* 2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq_file1%%.fastq*}_val_1_fastqc* 3_trimmed_fastqc
mv 2_trimgalore/${fq_file2%%.fastq*}_val_2_fastqc* 3_trimmed_fastqc

mkdir 0_rawfastq
mv $fq_file1 0_rawfastq
mv $fq_file2 0_rawfastq

## bismark alignment
mkdir 4_bismark_alignment
cd 4_bismark_alignment

## PE alignment
bismark --multicore 2 $genome_path -1 ../2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* -2 ../2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

## deduplicate BAMs (e.g. PCR duplicates)
deduplicate_bismark --bam -p ${fq_file1%%.fastq*}_val_1_bismark_bt2_pe.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# remove intermediate BAM
rm -v ${fq_file1%%.fastq*}_val_1_bismark_bt2_pe.bam

## PE methylation extraction
bismark_methylation_extractor --comprehensive --report --multicore 3 --buffer_size 8G \
	-p ${fq_file1%%.fastq*}_val_1_bismark_bt2_pe.deduplicated.bam \
       	2>&1 | tee -a ../${fileID}_logs_${dow}.log

#bedgraph creation on merged results
bismark2bedGraph --CX CpG*txt -o ${fileID}_CG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHG*txt -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHH*txt -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log


mkdir 5_output_files

mv 4_bismark_alignment/*.cov.gz 5_output_files

rm -v C*txt
rm -v *bed.gz

cd ../

mkdir 5_output_files
mv 4_bismark_alignment/*.cov.gz 5_output_files

echo "single C bed file with depth > 3"
## depth filter on single Cs
cd 5_output_files
zcat ${fileID}_CG.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

zcat ${fileID}_CHG.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

zcat ${fileID}_CHH.bed.gz.bismark.cov.gz | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7 = ($5 + $6)}' | awk '{ if ($7 >= 3) { print } }' > ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log

cd ../

echo "#####################"
echo "providing pipeline metrics to wgbs pipeline logfile..."
echo "#####################"

## get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/*_bismark_bt2_PE_report.txt | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/*_bismark_bt2_PE_report.txt | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/*_bismark_bt2_PE_report.txt | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/*_bismark_bt2_PE_report.txt | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

if [[ $fq_file1 == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file1 != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file1)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file1 == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4))
fi

if [[ $fq_file1 != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq*)
	flt_reads=$(($flt_reads / 4))
fi

#add it to the full pipeline logfile
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path}\t${type}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path}\t${type}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log
fi
