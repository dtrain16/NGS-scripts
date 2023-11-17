#!/usr/bin/env Rscript
# See BAM_to_EJC.sh to get input files (".5p.bed")
# Calculate 5'P reads counts 30-25 nt upstream of the 3' exon-exon junction
# Determine sum of 5P reads within EJC sites

library(tidyverse)

fls <- dir(pattern = ".5p.bed")

out <- NULL

for(i in fls){
        # Read in file
        input <- read.delim(i, head=F) %>%
        # Remove reads to plastid and mitochondria
        subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
        mutate(length = V7 - V6) %>%
        # features at least 50 bp in length
        subset(length > 49) %>%
        # calculate position relative to 3' end of feature
        mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3)) %>%
        subset(pos > -30 & pos < -25) %>%	
	group_by(V1, V2, V3, V8) %>% summarise(val=mean(V4)) %>%
	mutate(sample = sapply(strsplit(i, 'Aligned'), function(l) l[1]))
	
out <- rbind(input,out)

}

test <- group_by(out, sample, V8) %>%
	## sum read end depth across EJC window
	summarise(sum_val = sum(val)) %>%
	spread(sample, sum_val)

### NAs = 0 / no coverage
test[is.na(test)] <- 0

## Count table for NB-GLMs in edgeR
write.table(test, "ejc_counts.txt", sep='\t', quote=F, row.names=F)

