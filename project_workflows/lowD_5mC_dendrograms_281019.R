# making dendrograms for Eichten lowD brachy epigenomics
# libraries
library(tidyverse)
suppressPackageStartupMessages(library(dendextend))
library(pheatmap)

# read in metadata file
meta <- read.delim("lowD_metadata.txt") %>%
	select(SampleID, condition, PlantName, Accession, ClonalGroup)

my_grps <- data.frame(acc = unique(meta$Accession)) %>%
	mutate(Group = meta$ClonalGroup[match(acc, meta$Accession)]) %>%
	column_to_rownames("acc")

# read in tiled 5mC
i <- "CG_alltiles_merged_2017-05-05.txt"

## CG = "CG_alltiles_merged_2017-05-05.txt"
## CHG = "CHG_alltiles_merged.2017-05-04.txt"
# colnames(a)[4:ncol(a)] <- sapply(strsplit(colnames(a)[4:ncol(a)], ".fastq"), function(l) l[1])
# colnames(a)[4:ncol(a)] <- sapply(strsplit(colnames(a)[4:ncol(a)], ".S"), function(l) l[2])
# colnames(a)[4:ncol(a)] <- paste0("S",colnames(a)[4:ncol(a)])
## CHH = "CHH_alltiles_merged.2017-05-09.txt"

a <- read.delim(paste(i)) %>%
	gather(sample, met, -V1, -V2, -V3) %>%
	na.omit() %>%
	#mutate(condition = meta$condition[match(sample, meta$SampleID)]) %>%
	mutate(acc = meta$Accession[match(sample, meta$SampleID)]) %>%
	select(V1, V2, acc, met)

a <- group_by(a, V1, V2, acc) %>%
	summarise(avg_met = mean(met))

# memory allocation too great to pipe
a <- spread(a, acc, avg_met)

ann_colors <- list(
	Group = c(Bd21="coral", Clone1="royalblue", Clone2="darkgoldenrod1", 
	Clone3="darkolivegreen2", Clone4="darkorchid2", Clone5="forestgreen", 
	Clone6="firebrick1", Clone7="pink", HYB1="coral4", HYB2="bisque2")
)

# produce correlation matrix & heatmap
cor_matrix=as.matrix(cor(a[,3:ncol(a)],use='pairwise.complete.obs'))
pheatmap(cor_matrix, 
	cutree_cols = 3, 
	cutree_rows = 3,
	show_colnames = F,
	fontsize_row = 5,
	border_color = NA,
	annotation_colors = ann_colors,
	annotation_row = my_grps)

dev.off()

