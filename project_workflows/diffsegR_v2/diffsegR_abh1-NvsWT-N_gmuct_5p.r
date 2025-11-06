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
#conda install -n diffseg bioconda::r-scatterplot3d 
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
	sample    = c("abh1.N_1", "abh1.N_2", "abh1.N_3", "WT.N_1", "WT.N_2", "WT.N_3"),
	condition = c(rep("abh1.N", 3), rep( "WT.N", 3)),
	replicate = c(1:3,1:3),
	bam       = sapply(
		c("S15-5N_Aligned.sortedByCoord.out.bam", 
		"S9-20N_Aligned.sortedByCoord.out.bam",
		"S24-34N_Aligned.sortedByCoord.out.bam",
		"S5-3N_Aligned.sortedByCoord.out.bam",
		"S7-4N_Aligned.sortedByCoord.out.bam",
		"S11-10N_Aligned.sortedByCoord.out.bam"),
	function(bam) file.path(working_directory, bam)),
	isPairedEnd = rep(TRUE, 6),
	strandSpecific = rep(0, 6)	
)

#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_DERs <- NULL
out_segments <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## import data on experiment
data <- newExperiment(
	sampleInfo = sample_info,
	loci       = data.frame(seqid = i, chromStart = 1, chromEnd = stop),
	referenceCondition = "WT.N",
	otherCondition = "abh1.N",
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
	verbose = TRUE
)

#- extract DERs based on signifiance ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))
DERs <- subset(DERs, baseMean > 10)

out_DERs <- rbind(out_DERs,DERs)
}

## clear memory cache
gc()

out_DERs <- mutate(out_DERs, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, log2FoldChange, padj, baseMean)
write_tsv(out, "abh1-NvsWT-N_DERs_5p.bed", col_names=F)

