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
	sample    = c("WT.N_1","WT.N_2","WT.N_3","WT.N_4","WT.N_5","dxo1.N_1","dxo1.N_2","dxo1.N_3","dxo1.N_4","dxo1.N_5","dxo1.N_6","dxo1.N_7"),
	condition = c(rep("WT.N", 5), rep( "dxo1.N", 7)),
	replicate = c(1:5,1:7),
	bam = sapply( c("S5-3N_Aligned.sortedByCoord.out.bam","S7-4N_Aligned.sortedByCoord.out.bam", "S11-10N_Aligned.sortedByCoord.out.bam","S28-7N_Aligned.sortedByCoord.out.bam", "S30-8N_Aligned.sortedByCoord.out.bam", "S2-2N_Aligned.sortedByCoord.out.bam", "S17-6N_Aligned.sortedByCoord.out.bam", "S13-19N_Aligned.sortedByCoord.out.bam", "S36-45N_Aligned.sortedByCoord.out.bam","S37-46N_Aligned.sortedByCoord.out.bam","S38-47N_Aligned.sortedByCoord.out.bam","S39-48N_Aligned.sortedByCoord.out.bam"),
        function(bam) file.path(working_directory, bam)),
        isPairedEnd = rep(TRUE, 12),
        strandSpecific = rep(1, 12)
)


#- display sample information table -------------------------------------------#
knitr::kable(sample_info, row.names = FALSE)

genome <- read_tsv("~/ref_seqs/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.len", col_names=F)
genome <- subset(genome, X1 != "Mt" & X1 != "Pt")

### setup comparisons and loop for each chromosome
out_counts <- NULL

for(i in unique(genome$X1)){

chr <- paste(i)
stop <- genome$X2[genome$X1==i]

## import data on experiment
data <- newExperiment(
        sampleInfo 		= sample_info,
	referenceCondition 	= "WT.N",
        otherCondition  	= "dxo1.N",
        loci			= data.frame(seqid = i, chromStart = 1, chromEnd = stop),
	coverage 		= working_directory,        
	nbThreads  		= nb_threads,
        nbThreadsByLocus 	= nb_threads_locus
)

print(data)

## generate coverage profile from BAM
coverage(data = data, coverageType = "fivePrime", verbose = TRUE)

## transform coverage profile into per-base log2-FC and perform changepoint detection to define segments
features <- segmentationLFC(
	data  			= data, 
	modelSelectionType 	= "yao",
	alpha 			= 3,
	verbose 		= TRUE
)

## Quantify expression of segments
SExp <- counting(
	data 			 = data,
	features 		 = features,
	featureCountsType 	 = "fromBam",
	featureCountsOtherParams = list(read2pos = 5),
	verbose 		 = TRUE 
)

## subset to segments no longer than 10 nts
SExp_10 <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 11,]

## extract counts for segments
counts_10 <- SummarizedExperiment::assay(SExp_10)

out_counts <- rbind(out_counts,counts_10)

}

#clear memory cache
gc()

out_counts <- rownames_to_column(as.data.frame(out_counts), var="segment")
write_tsv(out_counts, "counts_full_Nsamples.tsv", col_names=T)

