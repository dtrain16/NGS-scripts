#!/usr/bin/env Rscript

## Installation
#### Conda
#conda create --name diffseg_v1
#conda install -n diffseg_v1 -c conda-forge r-tidyverse
#conda install -n diffseg_v1 -c bioconda bioconductor-deseq2
#conda install -n diffseg_v1 -c bioconda bioconductor-rsubread
#conda install -n diffseg_v1 -c bioconda bioconductor-rtracklayer
#conda install -n diffseg_v1 -c bioconda bioconductor-sparsematrixstats
#conda install -n diffseg_v1 -c bioconda bioconductor-delayedmatrixstats
#conda install -n diffseg_v1 -c conda-forge r-remotes
#conda activate diffseg_v1

#### R
#remotes::install_github("sanssouci-org/sanssouci")
#remotes::install_github("aLiehrmann/DiffSegR@f657435")

library(tidyverse)
library(DiffSegR)
nb_threads = 10
working_directory <- getwd()

#- create sample information table --------------------------------------------#
sample_info <- data.frame(
  sample    = c("abh1.N_1","abh1.N_2","abh1.N_3","abh1.C_1","abh1.C_2","abh1.C_3"),
  condition = c(rep("abh1.N", 3), rep( "abh1.C", 3)),
  replicate = c(1:3,1:3),
  bam       = sapply(
	c("S15-5N_Aligned.sortedByCoord.out.bam", "S9-20N_Aligned.sortedByCoord.out.bam", "S24-34N_Aligned.sortedByCoord.out.bam",
	"S16-5C_Aligned.sortedByCoord.out.bam", "S10-20C_Aligned.sortedByCoord.out.bam", "S25-34C_Aligned.sortedByCoord.out.bam"),
    function(bam) file.path(working_directory, bam)),
  coverage  = file.path(
    working_directory,
    paste0(c("abh1.N_1","abh1.N_2","abh1.N_3","abh1.C_1", "abh1.C_2", "abh1.C_3"), ".rds")))

#- save sample information table ----------------------------------------------# 
write.table(sample_info, file.path(working_directory, "sample_info.txt"))

#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_DERs <- NULL

for(i in unique(genome$X1)){
chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## load data
data <- loadData(
  sampleInfo   = file.path(working_directory,"sample_info.txt"),
  locus        = list(seqid = i, chromStart = 1, chromEnd = stop),
  referenceCondition = "abh1.C",
  isPairedEnd = TRUE,
  readLength = 150,
  coverageType = "fivePrime",
  stranded = FALSE,
  strandSpecific = 0,
  fromBam    = TRUE,
  nbThreads  = nb_threads,
  verbose = TRUE,
)

## Changepoint detection to define segments
SExp <- segmentation(
	data = data, 
	weightType = "unweighted", #zeroInflated : low counts have less weight
	modelSelectionType = "yao",
	featureCountsType = "fromBam",
	compressed = TRUE,
	alpha = 2,
	segmentNeighborhood = FALSE,
	Kmax = NULL,
	verbose = FALSE,
	nbThreadsGridSearch = 1,
	alphas = NULL,
	gridSearch = FALSE,
	outputDirectory = working_directory,
	nbThreadsFeatureCounts = nb_threads,
	strandSpecific = 0,
	read2pos = 5,
	isPairedEnd = TRUE
)

SExp_10 <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 11,]

dds <- dea(
	data      = data,
	SExp      = SExp_10,
	design    = ~condition,
	predicate = NULL,
	significanceLevel = 0.01,
	orderBy = "pvalue",
	verbose = TRUE
)

#- extract DERs based on signifiance ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$rejectedHypotheses,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))

out_DERs <- rbind(out_DERs,DERs)

}

#clear memory cache
gc()

out_DERs <- mutate(out_DERs, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, baseMean, baseVar, log2FoldChange, padj)
out <- subset(out, baseMean > 10)

write_tsv(out, "abh1-N_DERs_5p.bed", col_names=F)
