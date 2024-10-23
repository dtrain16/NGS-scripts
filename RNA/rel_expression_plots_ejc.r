#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_EJC.sh
# Computes normalized 5'P end frequency upstream of exon-exon junction (> 49 nt length, see Lee et al 2020 The Plant Cell)

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
	mutate(length = V7 - V6) %>%
# features at least 50 bp in length
	subset(length > 49) %>%
# calculate position relative to 3' end of feature
	mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3)) %>%
	mutate(rpm = V4 * rpm_scale)

# sum all reads in 50 nt window upstream of 3' end
exon_3p_sum <- subset(input, pos < 0 & pos > -51) %>%
	group_by(V8) %>%
	summarise(sum_rpm=sum(rpm))

# normalise depth per nt by sum of reads across 50 nt window and filter for raw read depth > 0
exon_3p <- subset(input, pos <= 0 & pos > -51) %>%
	mutate(sum_rpm = exon_3p_sum$sum_rpm[match(V8, exon_3p_sum$V8)]) %>%
	mutate(norm_counts = rpm/sum_rpm) %>%
	subset(abs(sum_rpm) > 1)

# Get sum of normalized reads (i.e. normalized occurrence of 5'P ends, Pi) then calculate relative frequency per nt
sum_exon_3p <- group_by(exon_3p, pos) %>%
	summarise(sum_norm_counts = sum(norm_counts), counts_raw = sum(V4), counts_rpm = sum(rpm)) %>% ## sum of normalized counts (Pi), raw coubts, and scaled counts (rpm) at each position
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

