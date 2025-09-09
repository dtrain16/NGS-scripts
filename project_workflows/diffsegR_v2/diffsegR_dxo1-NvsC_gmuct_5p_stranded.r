#!/usr/bin/env Rscript

# Installation
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
	sample    = c("dxo1.N_1","dxo1.N_2","dxo1.N_3","dxo1.N_4","WT.N_1","WT.N_2","WT.N_3","WT.N_4"),
	condition = c(rep("dxo1.N", 4), rep( "WT.N", 4)),
	replicate = c(1:4,1:4),
	bam = sapply( c("S5-3N_Aligned.sortedByCoord.out.bam","S7-4N_Aligned.sortedByCoord.out.bam", "S11-10N_Aligned.sortedByCoord.out.bam", "S30-8N_Aligned.sortedByCoord.out.bam",
"S36-45N_Aligned.sortedByCoord.out.bam","S37-46N_Aligned.sortedByCoord.out.bam","S38-47N_Aligned.sortedByCoord.out.bam","S39-48N_Aligned.sortedByCoord.out.bam"),
        function(bam) file.path(working_directory, bam)),
        isPairedEnd = rep(TRUE, 8),
        strandSpecific = rep(1, 8)
)


#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_DERs <- NULL
out_counts <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## import data on experiment
data <- newExperiment(
        sampleInfo = sample_info,
	referenceCondition = "WT.N",
        otherCondition  = "dxo1.N",
        loci		= data.frame(seqid = i, chromStart = 1, chromEnd = stop),
	coverage = working_directory,        
	nbThreads  = nb_threads,
        nbThreadsByLocus = nb_threads_locus
)

print(data)

## generate coverage profile from BAM
coverage(data = data, coverageType = "fivePrime", verbose = TRUE)

## transform coverage profile into per-base log2-FC and perform changepoint detection to define segments
features <- segmentationLFC(
	data  = data, 
	modelSelectionType = "yao",
	alpha = 2,
	verbose = TRUE
)

## subset to segments no longer than 10 nts
SExp_10 <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 11,]
## extract counts for segments
counts_10 <- SummarizedExperiment::assay(SExp_10)

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

out_counts <- rbind(out_counts,counts_10)
out_DERs <- rbind(out_DERs,DERs)

}

#clear memory cache
gc()


out_DERs <- mutate(out_DERs, derId = sapply(strsplit(featureId, "_"), function(l) paste0(l[1],":",l[2],"-",l[3])))
out <- select(out_DERs, seqnames, start, end, derId, baseMean, baseVar, log2FoldChange, padj)
out <- subset(out, baseMean > 10)

write_tsv(out, "dxo1vsWT_Nuc_DERs_5p.bed", col_names=F)

