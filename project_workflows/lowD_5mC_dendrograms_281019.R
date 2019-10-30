# making dendrograms for Eichten lowD brachy epigenomics

# libraries
library(tidyverse)
library(dendextend)
library(pheatmap)

# read in metadata file
meta <- read.delim("lowD_metadata.txt") %>%
	select(SampleID, condition, PlantName, Accession, ClonalGroup)

my_grps <- data.frame(acc = unique(meta$Accession)) %>%
	mutate(Group = meta$ClonalGroup[match(acc, meta$Accession)]) %>%
	column_to_rownames("acc")

# read in tiled 5mC
i <- "CG_alltiles_merged_2017-05-05.txt"

## CHG = "CHG_alltiles_merged.2017-05-04.txt"
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

# produce correlation matrix & heatmap
cor_matrix=as.matrix(cor(a[,3:ncol(a)],use='pairwise.complete.obs'))
pheatmap(cor_matrix, cutree_rows = 6, annotation = my_grps, fontsize = 4)
dev.off()

# hierarchical clustering and dendrogram
hc <- hclust(as.dist(1-cor(as.matrix(a[,3:ncol(a)]),use='pairwise.complete.obs')))
hc <- as.dendrogram(hc)
dendro_colors=as.numeric(meta$ClonalGroup)
labels_colors(hc) = dendro_colors[order.dendrogram(hc)]
plot(hc)
dev.off()

