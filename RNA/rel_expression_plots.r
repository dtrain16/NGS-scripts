#!/usr/bin/env Rscript
# args[1] = filename
# Run on output of BAM_to_wigs.sh

options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

# Read in file
input=read.delim(args[1],head=F)

# Remove plastids and unmatched rows
input=subset(input,input$V1!='ChrM' & input$V1!='ChrC' & input$V1 != 'Mt' & input$V1 != 'Pt')
input=subset(input,input[,ncol(input)] != -1)

# calculate normalized distance values for reads relative to feature
rel.dist=matrix(ifelse(input$V14==0,ifelse(input[,13]=="-",((input[,10] - (input[,2]))/(input[,10] - input[,9]))*1000,(((input[,2]) - input[,9])/(input[,10] - input[,9]))*1000),ifelse(input$V14>0,input$V14 + 1000,input$V14)),ncol=1)
input=cbind(input,rel.dist)
fixy=ifelse(input$rel.dist < 0 & input$V14==0,0,ifelse(input$rel.dist >1000 & input$V14==0,1000,input$rel.dist))
input$rel.dist=fixy

# bin read depth by distance
exp.bin=stats.bin(input$rel.dist,input$V4,N=100)
p.bin=cbind(matrix(exp.bin$centers,ncol=1),exp.bin$stats["mean",])
out=cbind(p.bin)
name <- sapply(strsplit(as.character(args[1]),'_'), function(l) l[1])
colnames(out)=c('pos',paste(name))
name2 <- sapply(strsplit(args[1], '\\.'), function(l) l[1])
write.table(out,paste(name2,'values.txt',sep=''),sep='\t', quote=F, row.names=F)
