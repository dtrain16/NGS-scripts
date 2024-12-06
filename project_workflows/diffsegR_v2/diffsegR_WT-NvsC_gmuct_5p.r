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
#conda install -n diffseg -c conda-forge r-remotes
#conda activate diffseg
#### R
#remotes::install_github("sanssouci-org/sanssouci")
#remotes::install_github("aLiehrmann/DiffSegR")

library(tidyverse)
library(DiffSegR)

## multi-threading options
nb_threads = 10
nb_threads_locus = 10

working_directory <- getwd()

#- create sample information table --------------------------------------------#
sample_info <- data.frame(
	sample    = c("WT.N_1", "WT.N_2", "WT.N_3", "WT.C_1", "WT.C_2", "WT.C_3"),
	condition = rep(c("WT.N", "WT.C"), each = 3),
	replicate = rep(1:3,2),
	bam       = sapply(
		c("S5-3N_Aligned.sortedByCoord.out.bam", 
		"S7-4N_Aligned.sortedByCoord.out.bam",
		"S11-10N_Aligned.sortedByCoord.out.bam",
		"S6-3C_Aligned.sortedByCoord.out.bam",
		"S8-4C_Aligned.sortedByCoord.out.bam",
		"S12-10C_Aligned.sortedByCoord.out.bam"
	), function(bam) file.path(working_directory, bam)),
	isPairedEnd = rep(TRUE, 6),
	strandSpecific = rep(0, 6)
)

#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

## genome file
genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_DERs <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## import data on experiment
data <- newExperiment(
	sampleInfo = sample_info,
	loci       = data.frame(seqid = i, chromStart = 1, chromEnd = stop),
	referenceCondition = "WT.C",
	otherCondition = "WT.N",
	nbThreads  = nb_threads,
	nbThreadsByLocus = nb_threads_locus,
	coverage = working_directory
)

print(data)

## generate coverage profile from BAM
coverage(data = data, coverageType = "fivePrime", verbose = TRUE)

## transform coverage profile into per-base log2-FC and perform changepoint detection to define segments
features <- segmentationLFC(
	data  = data, 
	alpha = 2,
	modelSelectionType = "yao",
	verbose = TRUE
)

## Quantify expression of segments
SExp <- counting(
	data = data,
	features = features,
	featureCountsType = "fromBam",
	featureCountsOtherParams = list(read2pos = 5),
	verbose = TRUE 
)

#- subset to segments with width < 11 nt --------------------------------------#
SExp_10 <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 11,]


# differential exprssion analysis
dds <- dea(
	SExp        = SExp_10, 
	design      = ~condition,
	significanceLevel = 0.01,
	verbose = TRUE,
	predicate = NULL	
)

#- extract DERs based on signifiance ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))

out_DERs <- rbind(out_DERs,DERs)
}

## clear memory cache
gc()

out_DERs <- mutate(out_DERs, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, baseMean, baseVar, maxCooks, log2FoldChange, padj)
out$cov <- sqrt(out$baseVar)/out$baseMean
out <- subset(out, baseMean > 10 & cov < 1)

write_tsv(out, "WT-N_DERs_5p.bed", col_names=F)

