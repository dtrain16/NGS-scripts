#!/usr/bin/env Rscript
# See BAM_to_STOP.sh to get input files ("50bp.5p.bed")
# Determine transripts with differential co-translational decay signal, inferred by abundance of 5'P reads 17-20 nt upstream of the STOP codon.

library(tidyverse)

fls <- dir(pattern = "_50bp.5p.bed")

out <- NULL

for(i in fls){
        # Read in file
        input <- read.delim(i, head=F) %>%
        # Remove reads to plastid and mitochondria
        subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
        # calculate position relative to 3' end of feature
        mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3)) %>%
        subset(pos > -20 & pos < -14) %>%	
	group_by(V1, V2, V3, V9) %>% summarise(val=mean(V4)) %>%
	mutate(sample = sapply(strsplit(i, 'Aligned'), function(l) l[1]))
	
out <- rbind(input,out)

}

test <- group_by(out, sample, V9) %>%
	## sum read end depth across CTRD window
	summarise(sum_val = sum(val)) %>%
	spread(sample, sum_val)

### NAs = 0 / no coverage
test[is.na(test)] <- 0

## Count table for NB-GLMs in edgeR
write.table(test, "ctrd_rna_dxo1_counts.txt", sep='\t', quote=F, row.names=F)


