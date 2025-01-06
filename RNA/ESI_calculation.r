#!/usr/bin/env Rscript
# args[1] = filename
# Required for BAM_to_ESI.sh
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
        mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3)) %>%
	 subset(pos > -50 & pos < 1)

# get total end counts in window
stop_5p_sum <- group_by(input, V8) %>%
	summarise(avg_frame=mean(V4), avg_frame_true=sum(V4)/50, total_counts=sum(V4))

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
a1 <- group_by(input, V8) %>% 
	subset(pos == -28 | pos == -27) %>%
	summarise(avg_ejc = mean(V4), avg_ejc_true = sum(V4)/2) %>%
	mutate(avg_frame = stop_5p_sum$avg_frame[match(V8, stop_5p_sum$V8)]) %>%
	mutate(total_counts = stop_5p_sum$total_counts[match(V8, stop_5p_sum$V8)]) %>%
	subset(total_counts >= 10) %>%
	mutate(esi = avg_ejc/avg_frame)

# name output filea
name <- sapply(strsplit(as.character(args[1]),'Aligned'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,"ESI.pdf"))
plot(y=log2(a1$avg_ejc), x=log2(a1$avg_frame), col = ifelse(a1$tsi > 2, "salmon", "grey"))
abline(a=0, b=1)
dev.off()

## output
write.table(a1, paste0(name,"ESI.txt"), sep='\t', quote=F, row.names=F)

