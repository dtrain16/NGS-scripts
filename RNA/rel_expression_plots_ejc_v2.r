#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_EJC_3p.sh
# Computes normalized 5'P end frequency for last 50 nt upstream of exon-exon junction (> 49 nt length, see Lee et al 2019 Plant Cell)

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input=read.delim(args[1],head=F) %>% 
# Remove plastids and unmatched rows
	subset(input$V1!='ChrM' & input$V1!='ChrC' & input$V1 != 'Mt' & input$V1 != 'Pt') %>%
	mutate(length = V7 - V6) %>%
# features at least 50 bp in length
	subset(length > 49) %>%
# calculate position relative to 3' end of feature
	mutate(pos_3p = ifelse(V10 == "+", V2-V7, V6-V3))


# sum all reads in 50 nt window upstream of 3' end
exon_3p_sum <- subset(input, pos_3p < 0 & pos_3p > -50) %>%
	group_by(V8) %>%
	summarise(reads=sum(V4))

# normalise depth per nt by sum of reads across 50 nt window
exon_3p <- subset(input, pos_3p < 0 & pos_3p > -51) %>%
	mutate(sum_reads = exon_3p_sum$reads[match(V8, exon_3p_sum$V8)]) %>%
	mutate(norm_reads = V4/sum_reads) %>%
	subset(sum_reads > 0)

# Get sum of normalized reads (i.e.normalized occurrence of 5'P ends [Pi] in Lee et al 2019 Plant Cell) then calculate relative frequency per nt
sum_exon_3p <- group_by(exon_3p, pos_3p) %>% summarise(sum_norm_reads = sum(norm_reads)) %>% mutate(total_reads = sum(sum_reads)) %>% mutate(rel_freq = sum_norm_reads/total_reads)

name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".pdf"))
plot(y=sum_exon_3p$rel_freq, x= sum_exon_3p$pos_3p)
dev.off()

## output
write.table(sum_exon_3p, paste0(name,".5p.txt"), sep='\t', quote=F, row.names=F)

