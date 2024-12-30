#!/usr/bin/env Rscript
# Run on output of graft_nad_adprc_enrich.sh

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input <- read.delim(args[1],head=F) %>% 
# Remove reads to plastid and mitochondria
	subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
# calculate position relative to 3' end of feature
	mutate(pos = ifelse(V10 == "+", V7-V2, V3-V6))

# normalise depth per nt by sum of reads across 50 nt window and filter for raw read depth > 0
exon_3p <- subset(input, pos <= 0 & pos > -51) %>%
	mutate(sum_rpm = exon_3p_sum$sum_rpm[match(V8, exon_3p_sum$V8)]) %>%
	mutate(norm_counts = rpm/sum_rpm) %>%
	subset(abs(sum_rpm) > 1)

# Get sum of normalized reads (i.e. normalized occurrence of 5'P ends, Pi) then calculate relative frequency per nt
sum_exon_3p <- group_by(exon_3p, pos) %>%
	summarise(sum_norm_counts = sum(norm_counts), counts_raw = sum(V4), counts_rpm = sum(rpm)) %>% ## sum of normalized counts (Pi), raw coubts, and scaled counts (rpm) per nt
	mutate(total_counts = sum(sum_norm_counts)) %>% ## sum of all normalized counts 
	mutate(rel_freq = sum_norm_counts/total_counts) %>% ## freq of normalized counts per position relative to all normalized counts
	select(pos, sum_norm_counts, rel_freq) # take columns for Pi and relative frequency of Pi

# name output file
name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".3p.pdf"))
plot(y=sum_exon_3p$rel_freq, x= sum_exon_3p$pos)
dev.off()

## output
write.table(sum_exon_3p, paste0(name,".ejc.txt"), sep='\t', quote=F, row.names=F)

