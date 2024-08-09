#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_STOP.sh
# Calculate EJC stalling index from 5'P end counts upstream of exon-exon junctions

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input <- read.delim(args[1],head=F) %>% 
# Remove reads to plastid and mitochondria
	subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
#exon at least 50 nucleotides
	mutate(length = V7 - V6) %>%
        subset(length > 49) %>%
# calculate position relative to 3' end of feature
        mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3))

# get total end counts in window
stop_5p_sum <- group_by(input, V8) %>%
	summarise(avg_frame=mean(V4), avg_frame_true=sum(V4)/100, total_counts=sum(V4))

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
a1 <- group_by(input, V8) %>% 
	subset(pos == -16 | pos == -17) %>%
	summarise(avg_ctrd = mean(V4), avg_ctrd_true = sum(V4)/2) %>%
	mutate(avg_frame = stop_5p_sum$avg_frame[match(V8, stop_5p_sum$V8)]) %>%
	mutate(total_counts = stop_5p_sum$total_counts[match(V8, stop_5p_sum$V8)]) %>%
	subset(total_counts >= 20) %>%
	mutate(tsi = avg_ctrd/avg_frame)

# name output filea
name <- sapply(strsplit(as.character(args[1]),'Aligned'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,"_TSI.pdf"))
plot(y=log2(a1$avg_ctrd), x=log2(a1$avg_frame))
abline(a=0, b=1)
dev.off()

## output
write.table(a1, paste0(name,"_TSI.txt"), sep='\t', quote=F, row.names=F)

