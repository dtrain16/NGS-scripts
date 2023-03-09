#!/usr/bin/env Rscript
# See BAM_to_EJC.sh to get input files ("10bp.50.bed")
# Determine transripts with differential occupancy of EJC complex, inferred by abundance of 5'P reads 30-25 nt upstream of the 3' end of an exon.

library(fields)
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

test <- mutate(out, genotype = sapply(strsplit(sample, "_"), function(l) l[1])) %>%
	group_by(sample, genotype, V8) %>%
	## sum read end depth
	summarise(sum_val = sum(val)) %>%
	spread(sample, sum_val)

### NAs = 0 / no coverage
test[is.na(test)] <- 0
## at least 3 read ends in 2 samples
keep <- rowSums(test[3:6] > 3) > 1
y <- test[keep,]

y1 <- 	gather(y,sample, val, -genotype, -V8) %>%
	group_by(genotype, V8) %>%
	# mean read ends per genotype
	summarise(mean_val = mean(val)) %>% 
	spread(genotype, mean_val)

test[is.na(test)] <- 0
test <- mutate(test, fc = (dxo1+1)/(WT+1))
median(test$fc)
ejc_out <- subset(test, fc > median(test$fc)) 

## output
write.table(ejc_out, "ejc_rna_dxo1.txt", sep='\t', quote=F, row.names=F)

