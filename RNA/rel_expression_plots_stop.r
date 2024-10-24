#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_STOP.sh
# computes 5'P end depth relative to first nucleotide of STOP or START codon 

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

rpm_scale <- as.numeric(paste(args[2]))

# Read in file
input <- read.delim(args[1],head=F) %>% 
# Remove reads to plastid and mitochondria
	subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
# calculate position relative to first base of STOP or START codon
	mutate(pos = ifelse(V10 == "+", V2-V6, V7-V3)) %>%
	mutate(rpm = V4 * rpm_scale)

# sum all reads in 50 nt window upstream of 3' end
stop_5p_sum <- group_by(input, V8) %>%
	summarise(sum_rpm=sum(rpm), sum_reads=sum(V4))

# normalise depth per nt by sum of reads across window
stop_5p <- mutate(input, sum_rpm = stop_5p_sum$sum_rpm[match(V8, stop_5p_sum$V8)]) %>%
	mutate(norm_counts = rpm/sum_rpm) %>%
	subset(sum_rpm > 1)

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends, Pi) then calculate relative frequency per nt
sum_stop_5p <- group_by(stop_5p, pos) %>% 
	summarise(sum_norm_counts = sum(norm_counts), counts_raw = sum(V4), counts_rpm = sum(rpm)) %>% ## sum of normalized counts (Pi), raw counts, and scaled counts (rpm) per nt  
	mutate(total_counts = sum(sum_norm_counts)) %>% ## sum of all normalized counts 
	mutate(rel_freq = sum_norm_counts/total_counts) %>% ## freq of normalized counts per position relative to all normalized counts
	select(pos, sum_norm_counts, rel_freq)

# name output file
name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".stop.pdf"))
plot(y=sum_stop_5p$rel_freq, x= sum_stop_5p$pos)
dev.off()

## output
write.table(sum_stop_5p, paste0(name,".stop.txt"), sep='\t', quote=F, row.names=F)

