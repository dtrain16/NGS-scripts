#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_EJC.sh
# Computes normalized 5'P end frequency across 50 nt at 5' of exon (> 49 nt length, see Lee et al 2019 Plant Cell)

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input <- read.delim(args[1],head=F) %>% 
# Remove reads to plastid and mitochondria
	subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
	mutate(length = V7 - V6) %>%
# features at least 50 bp in length
	subset(length > 49) %>%
# calculate position relative to 3' end of feature
	mutate(pos_5p = ifelse(V10 == "+", V2-V6, V7-V3))

# sum all reads in 50 nt window upstream of 3' end
exon_5p_sum <- subset(input, pos_5p > 0 & pos_5p < 51) %>%
	group_by(V8) %>%
	summarise(reads=sum(V4))

# normalise depth per nt by sum of reads across 50 nt window
exon_5p <- subset(input, pos_5p > 0 & pos_5p < 51) %>%
	mutate(sum_reads = exon_5p_sum$reads[match(V8, exon_5p_sum$V8)]) %>%
	mutate(norm_reads = V4/sum_reads) %>%
	subset(abs(sum_reads) > 0)

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
sum_exon_5p <- group_by(exon_5p, pos_5p) %>% summarise(sum_norm_reads = sum(norm_reads)) %>% mutate(total_reads = sum(sum_norm_reads)) %>% mutate(rel_freq = sum_norm_reads/total_reads)

# name output file
name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".5p.pdf"))
plot(y=sum_exon_5p$rel_freq, x= sum_exon_5p$pos_5p)
dev.off()

## output
write.table(sum_exon_3p, paste0(name,".5p.txt"), sep='\t', quote=F, row.names=F)

