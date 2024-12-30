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

## multi-threading options
nb_threads = 10

working_directory <- getwd()

#- create sample information table --------------------------------------------#
sample_info_bdg <- data.frame(
	sample    = c("WT.N_1", "WT.N_2", "WT.N_3", "BDG_1", "BDG_2"),
	condition = c(rep("WT.N", each = 3),rep("BDG", each=2)),
	replicate = c(1:3,1:2),
	bam       = sapply(
		c("S5-3N_R1_Aligned.sortedByCoord.out.bam", "S7-4N_R1_Aligned.sortedByCoord.out.bam", "S11-10N_R1_Aligned.sortedByCoord.out.bam",
		"WT-BDG_rep1_Aligned.sortedByCoord.out.bam", "WT-BDG_rep2_Aligned.sortedByCoord.out.bam"),
	function(bam) file.path(working_directory, bam)),
	coverage  = file.path(working_directory, paste0(c("WT.N_1", "WT.N_2", "WT.N_3", "BDG_1", "BDG_2"), ".rds")))

sample_info_hmc  <- data.frame(
        sample    = c("WT.N_1", "WT.N_2", "WT.N_3", "HMC_1", "HMC_2"),
        condition = c(rep("WT.N", each = 3),rep("HMC",each=2)),
        replicate = c(1:3,1:2),
        bam       = sapply(
                c("S5-3N_R1_Aligned.sortedByCoord.out.bam", "S7-4N_R1_Aligned.sortedByCoord.out.bam", "S11-10N_R1_Aligned.sortedByCoord.out.bam",
                "WT-HMC_rep1_Aligned.sortedByCoord.out.bam", "WT-HMC_rep2_Aligned.sortedByCoord.out.bam"),
        function(bam) file.path(working_directory, bam)),
        coverage  = file.path(working_directory, paste0(c("WT.N_1", "WT.N_2", "WT.N_3", "HMC_1", "HMC_2"), ".rds")))

#- save sample information table ----------------------------------------------# 
write.table(sample_info_bdg, file.path(working_directory, "sample_info_bdg.txt"))

write.table(sample_info_hmc, file.path(working_directory, "sample_info_hmc.txt"))

#- display sample information table -------------------------------------------#
knitr::kable(sample_info_bdg, row.names = FALSE)
knitr::kable(sample_info_hmc, row.names = FALSE)

## genome file
genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_DERs_bdg <- NULL
out_DERs_hmc <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## import data on experiment
data1 <- loadData(sampleInfo   = file.path(working_directory,"sample_info_bdg.txt"),
	locus        = list(seqid = i, chromStart = 1, chromEnd = stop),
	referenceCondition = "BDG",
	isPairedEnd = FALSE,
	readLength = 50,
	coverageType = "threePrime",
	stranded = FALSE,
	strandSpecific = 0,
	fromBam    = TRUE,
	nbThreads  = nb_threads,
	verbose = TRUE,
)

data2 <- loadData(sampleInfo   = file.path(working_directory,"sample_info_hmc.txt"),
        locus      = list(seqid = i, chromStart = 1, chromEnd = stop),
        referenceCondition = "HMC",
        isPairedEnd = FALSE,
        readLength = 50,
        coverageType = "threePrime",
        stranded = FALSE,
        strandSpecific = 0,
        fromBam    = TRUE,
        nbThreads  = nb_threads,
        verbose = TRUE,
)


## Changepoint detection to define segments
SExp1 <- segmentation(data = data1, 
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
	read2pos = 3,
	isPairedEnd = FALSE
)

SExp2 <- segmentation(data = data2,
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
        read2pos = 3,
        isPairedEnd = FALSE
)


#- subset to segments with width < 11 nt --------------------------------------#
SExp1_10 <- SExp1[as.data.frame(SummarizedExperiment::rowRanges(SExp1))$width < 11,]
SExp2_10 <- SExp2[as.data.frame(SummarizedExperiment::rowRanges(SExp2))$width < 11,]

# differential exprssion analysis
dds <- dea(data = data1,
	SExp    = SExp1_10,
	design  = ~condition,
	predicate = NULL,
	significanceLevel = 0.01,
	verbose = TRUE
)

#- extract DERs based on signifiance ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))

out_DERs_bdg <- rbind(out_DERs_bdg, DERs)

# differential exprssion analysis
dds <- dea(data = data2,
        SExp    = SExp2_10,
        design  = ~condition,
        predicate = NULL,
        significanceLevel = 0.01,
        verbose = TRUE
)

#- extract DERs based on signifiance ----------------------------------------#
DERs <- dds[SummarizedExperiment::mcols(dds)$rejectedHypotheses,]
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))

out_DERs_hmc <- rbind(out_DERs_hmc, DERs)

}

## clear memory cache
gc()

out_DERs <- mutate(out_DERs_bdg, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, baseMean, baseVar, log2FoldChange, padj)
out <- subset(out, baseMean > 10)

write_tsv(out, "WT-NvsBDG_5p.bed", col_names=F)

out_DERs <- mutate(out_DERs_hmc, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, baseMean, baseVar, log2FoldChange, padj)
out <- subset(out, baseMean > 10)

write_tsv(out, "WT-NvsHMC_5p.bed", col_names=F)








