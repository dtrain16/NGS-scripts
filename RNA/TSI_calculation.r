#!/usr/bin/env Rscript
# args[1] = filename
# Runs with BAM_to_TSI.sh
# Calculate terminal stalling index from 5'P end counts surrounding the stop codon

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input <- read.delim(args[1],head=F) %>% 
# Remove reads to plastid and mitochondria
	subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
# calculate position relative to first base of stop codon
	mutate(pos = ifelse(V10 == "+", V2-V6, V7-V3)) %>%
	subset(pos > -50 & pos < 1)

# get total end counts in window
stop_5p_sum <- group_by(input, V8) %>%
	summarise(avg_frame=mean(V4), avg_frame_true=sum(V4)/50, total_counts=sum(V4))

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
a1 <- group_by(input, V8) %>% 
	subset(pos == -16 | pos == -17) %>%
	summarise(avg_ctrd = mean(V4), avg_ctrd_true = sum(V4)/2) %>%
	mutate(avg_frame = stop_5p_sum$avg_frame[match(V8, stop_5p_sum$V8)]) %>%
	mutate(total_counts = stop_5p_sum$total_counts[match(V8, stop_5p_sum$V8)]) %>%
	subset(total_counts >= 10) %>%
	mutate(tsi = avg_ctrd/avg_frame)

# name output filea
name <- sapply(strsplit(as.character(args[1]),'Aligned'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,"TSI.pdf"))
plot(y=log2(a1$avg_ctrd), x=log2(a1$avg_frame), col = ifelse(a1$tsi > 2, "salmon", "grey"))
abline(a=0, b=1)
dev.off()

## output
write.table(a1, paste0(name,"TSI.txt"), sep='\t', quote=F, row.names=F)

