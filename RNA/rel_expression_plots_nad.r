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
# calculate position relative to 5' end of feature
	mutate(pos = ifelse(V12 == "+", V2-V8, V9-V3)) %>%
	subset(pos <= 50 & pos >= -50) %>%
	subset(V4 >= 1) ## at least 1 RPM in ADPRC sample

# summarise prop NAD+
df1 <- group_by(input, pos) %>%
	summarise(mean_prop = mean(V6)) # mean prop NAD = A+ / A-

# name output file
name <- sapply(strsplit(as.character(args[1]),'.bed'), function(l) l[1])

## diagnostic plot on single sample
pdf(paste0(name,".pdf"))
plot(y=df1$mean_prop, x= df1$pos)
dev.off()

## output
write.table(df1, paste0(name,".ejc.txt"), sep='\t', quote=F, row.names=F)

