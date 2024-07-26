#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_STOP.sh
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
	subset(pos > -101 & pos < 1)

# sum all reads in 50 nt window upstream of 3' end
stop_5p_sum <- group_by(input, V8) %>%
	summarise(reads=sum(V4))

# normalise depth per nt by sum of reads across window
stop_5p <- mutate(input, sum_reads = stop_5p_sum$reads[match(V8, stop_5p_sum$V8)]) %>%
	subset(abs(sum_reads) > 9)

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
a1 <- group_by(stop_5p, V8) %>% 
	subset(pos == -16 | pos == -17) %>%
	summarise(raw_counts = mean(V4)) 

# name output filea
name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".stop.pdf"))
plot(y=sum_stop_5p$rel_freq, x= sum_stop_5p$pos)
dev.off()

## output
write.table(sum_stop_5p, paste0(name,".stop.txt"), sep='\t', quote=F, row.names=F)

