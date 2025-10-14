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
	sample    = c("WT.A_1","WT.A_2", "N.A_1","N.A_2","WT.X_1","WT.X_2","N.X_1","N.X_2"),
	condition = c(rep("ADPRC+", 4), rep( "ADPRC-", 4)),
	replicate = c(1:4,1:4),
	bam = sapply( c("WT_A_rep1_Aligned.sortedByCoord.out.bam", "WT_A_rep2_Aligned.sortedByCoord.out.bam", "N_A_rep1_Aligned.sortedByCoord.out.bam", "N_A_rep2_Aligned.sortedByCoord.out.bam",
"WT_X_rep1_Aligned.sortedByCoord.out.bam", "WT_X_rep2_Aligned.sortedByCoord.out.bam","N_X_rep1_Aligned.sortedByCoord.out.bam","N_X_rep2_Aligned.sortedByCoord.out.bam"),
        function(bam) file.path(working_directory, bam)),
        isPairedEnd = rep(FALSE, 8),
        strandSpecific = rep(1, 8)
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
	referenceCondition 	= "ADPRC+",
        otherCondition  	= "ADPRC-",
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
SExp_filter <- SExp[as.data.frame(SummarizedExperiment::rowRanges(SExp))$width < 5,]

## extract counts for segments
counts <- SummarizedExperiment::assay(SExp_filter)

out_counts <- rbind(out_counts,counts)

}

#clear memory cache
gc()

out_counts <- rownames_to_column(as.data.frame(out_counts), var="segment")
write_tsv(out_counts, "Graft-NAD_WTvsnudt19_counts.tsv", col_names=T)

