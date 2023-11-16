#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_bedgraph.sh or BAM_to_bedgraph_5p.sh
# calculate total read depth across features of interest

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input <- read.delim(args[1],head=F)

input <- subset(input,input[,ncol(input)] != -1)

# calculate length of feature and per million scaling factor (total reads)
input <- subset(input, V1 !='ChrM' & V1!='ChrC' & V1 != 'Mt' & V1 != 'Pt') %>% # Remove plastids and unmatched rows
	mutate(input, length = V7 - V6)
	
# Determine total read depth per feature then calculate RPM using per million scaling factor
out <- group_by(input, V8, length) %>%
	summarise(read_depth = sum(V4))

name <- sapply(strsplit(args[1], 'bed'), function(l) l[1])
write.table(out, paste(name,'feature_depth.txt',sep=''), sep='\t', quote=F, row.names=F)

