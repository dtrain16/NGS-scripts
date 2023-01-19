#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_wigs.sh
# Computes normalized 5'P read frequency across full exon length

options(echo=T)
library(fields)
library(tidyverse)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input=read.delim(args[1],head=F)

# Remove plastids and unmatched rows
input <- subset(input,V1!='ChrM' & V1!='ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
	subset(input[,ncol(input)] != -1) %>%
	mutate(length = V7 - V6) %>%
        subset(length > 49) 

# calculate normalize 5'P occurence based on total reads per exon (only within the exon)
exon_reads <- group_by(input, V8) %>%
	summarise(reads=sum(V4))

input <- mutate(input, exon_reads = exon_reads$reads[match(V8, exon_reads$V8)]) %>%
	subset(exon_reads != 0) %>%
	mutate(norm_reads = V4/exon_reads) %>%
	mutate(norm_reads = ifelse(exon_reads < 0, norm_reads * -1, norm_reads))

# calculate normalized distance values for reads relative to feature
rel.dist=matrix(ifelse(input$V10=="-",(input$V2 - input$V6),(input$V3 - input$V7)),ncol=1)

rel.dist=matrix(ifelse(input$V11==0,ifelse(input[,10]=="-",((input[,7] - (input[,2]))/(input[,7] - input[,6]))*100,(((input[,2]) - input[,6])/(input[,7] - input[,6]))*100),ifelse(input$V11>0,input$V11 + 100,input$V11)),ncol=1)

input=cbind(input,rel.dist)

fixy=ifelse(input$rel.dist < 0 & input$V11==0,0,ifelse(input$rel.dist > 100 & input$V11==0 , 100, input$rel.dist))
input$rel.dist=fixy

# bin read depth by distance
exp.bin=stats.bin(input$rel.dist,input$norm_reads,N=100)
p.bin=cbind(matrix(exp.bin$centers,ncol=1),exp.bin$stats["mean",])
out=cbind(p.bin)
name <- sapply(strsplit(as.character(args[1]),'_'), function(l) l[1])
colnames(out)=c('pos',paste(name))
name2 <- sapply(strsplit(args[1], '\\.'), function(l) l[1])
name3 <- sapply(strsplit(args[1], '_'), function(l) l[3])
name3 <- gsub(".bed",".values.txt",name3)
write.table(out,paste(name2,name3,sep='.'),sep='\t', quote=F, row.names=F)


