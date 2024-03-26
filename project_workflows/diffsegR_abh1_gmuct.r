#!/usr/bin/env Rscript

## Installation
#### Conda
#conda create --name diffseg
#conda install -n diffseg -c conda-forge r-tidyverse
#conda install -n diffseg -c bioconda bioconductor-deseq2
#conda install -n diffseg -c bioconda bioconductor-rsubread
#conda install -n diffseg -c bioconda bioconductor-rtracklayer
#conda install -n diffseg -c bioconda bioconductor-sparsematrixstats
#conda install -n diffseg -c bioconda bioconductor-delayedmatrixstats
#conda activate diffseg
#### R
#remotes::install_github("sanssouci-org/sanssouci")
#remotes::install_github("aLiehrmann/DiffSegR")

library(tidyverse)
library(DiffSegR)
nb_threads = 8

working_directory <- getwd()

#- create sample information table --------------------------------------------#
sample_info <- data.frame(
  sample    = c("abh1_1", "abh1_2", "WT_1", "WT_2"),
  condition = rep(c("dxo1", "WT"), each = 2),
  replicate = rep(1:2,2),
  bam       = sapply(
    c("abh1_rep1Aligned.sortedByCoord.out.bam", 
      "abh1_rep2Aligned.sortedByCoord.out.bam",
      "WT_rep1Aligned.sortedByCoord.out.bam",
      "WT_rep2Aligned.sortedByCoord.out.bam"
    ),
    function(bam) file.path(working_directory, bam)
  ),
  coverage  = file.path(
    working_directory,
    paste0(c("abh1_1", "abh1_2", "WT_1", "WT_2"), ".rds")
  )
)

#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### WT vs abh1
## repeat DER identification for each chromosome

out_DERs <- NULL
out_segments <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## load data
data <- loadData(
  sampleInfo   = file.path(working_directory,"sample_info.txt"),
  locus        = list(seqid = i, chromStart = 1, chromEnd = stop),
  referenceCondition = "WT",
  isPairedEnd = FALSE,
  readLength = 50,
  coverageType = "fivePrime",
  stranded = FALSE,
  #strandSpecific = 0,
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
	read2pos = 5, #NULL, 5, 3
	isPairedEnd = FALSE
)

#- subset to segments with width < 5 nt ------------------------------------------------#
segment_coordinates <- as.data.frame(SummarizedExperiment::rowRanges(SExp))
segment_coordinates <- subset(segment_coordinates, width < 6)
SExp_5 <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 6,]

# differential exprssion analysis
dds <- dea(
  data              = data,
  SExp              = SExp_5, ##on segments < 5 nt width
  design            = ~condition,
  sizeFactors       = NA,
  significanceLevel = 0.05,
  predicate = NULL,
  postHoc_significanceLevel = 0.05,
  postHoc_tdpLowerBound = 0.95,
  orderBy = "pvalue",
  verbose = FALSE,
  dichotomicSearch = FALSE
)

#- display DERs by position ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$rejectedHypotheses,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))
dim(DERs) #number

out_segments <- rbind(out_segments, segment_coordinates)
out_DERs <- rbind(out_DERs,DERs)
}

## clear memory cache
gc()

pdf("qc_plots_abh1.pdf")
hist(out_DERs$modelMean)
hist(out_DERs$log2FoldChange)
plot(x=out_DERs$modelMean, y=out_DERs$log2FoldChange)
plot(x=out_DERs$log2RefMean, y=out_DERs$log2OtherMean)
dev.off()

out2 <- select(out_DERs, featureId, seqnames, start, end, width, strand, log2FoldChange, padj, log2RefMean, log2OtherMean)
write_tsv(out2, "abh1_DERs.tsv", col_names=F)

out3 <- select(out2, seqnames, start, end, featureId, log2FoldChange, strand)
write_tsv(out3, "abh1_DERs.bed", col_names=F)

knitr::kable(head(out_DERs[colnames(out_DERs)%in%c("seqnames","start","end","width","strand","log2FoldChange","padj","log2RefMean","log2OtherMean")]))


