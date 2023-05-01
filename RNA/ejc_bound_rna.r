#!/usr/bin/env Rscript
# See BAM_to_EJC.sh to get input files ("10bp.50.bed")
# Determine transripts with differential occupancy of EJC complex, inferred by abundance of 5'P reads 30-25 nt upstream of the 3' end of an exon.

library(tidyverse)

fls <- dir(pattern = "_10bp.5p.bed")

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
        subset(pos > -31 & pos < -24) %>%	
	group_by(V1, V2, V3, V8) %>% summarise(val=mean(V4)) %>%
	mutate(sample = sapply(strsplit(i, 'Aligned'), function(l) l[1]))
	
out <- rbind(input,out)

}

test <- group_by(out, sample, V8) %>%
	## sum read end depth
	summarise(sum_val = sum(val)) %>%
	spread(sample, sum_val)

### NAs = 0 / no coverage
test[is.na(test)] <- 0
## at least 3 read ends in 2 samples
keep <- rowSums(test[2:5] > 5) > 1
y <- test[keep,]

#y1 <- 	gather(y, sample, val, -V8) %>%
#	mutate(genotype = sapply(strsplit(sample, "_"), function(l) l[1])) %>%
#	group_by(genotype, V8) %>%
	# mean read ends per genotype
#	summarise(mean_val = mean(val)) %>% 
#	spread(genotype, mean_val) %>%
#	mutate(fc = (dxo1+1)/(WT+1))

## output
write.table(y, "ejc_rna_dxo1.txt", sep='\t', quote=F, row.names=F)


