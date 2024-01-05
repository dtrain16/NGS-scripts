#!/usr/bin/env Rscript
# See BAM_to_bedgraph_5p.sh to get input files (".5p.bed")
# Determine abundance of max 5P peak (Max 5P) across features of interest

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
	group_by(V1, V2, V3, V8) %>% 
	summarise(val=mean(V4)) %>%
	mutate(sample = sapply(strsplit(i, 'Aligned'), function(l) l[1]))
	
out <- rbind(input,out)

}

test <- group_by(out, sample, V8) %>%
	## get maximum 5'P end depth within EJC window
	summarise(max5p = max(val)) %>%
	spread(sample, max5p)

### NAs = 0 / no coverage
test[is.na(test)] <- 0

## Count table for NB-GLMs in edgeR
write.table(test, "max5p.txt", sep='\t', quote=F, row.names=F)

