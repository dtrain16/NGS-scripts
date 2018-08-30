#!/usr/bin/env Rscript
# Produce methylome correlation matrix using specified BED files

args = commandArgs(trailingOnly=T)

if(length(args) != 4){
	print('Args missing!')
	print('USAGE: merge_wigs.r <context> <bin size> <cov> <type [1,reps vs 2,sum]>')
	print('EXAMPLE merge_wigs.r CHH 0 0 1')
	quit()
	n
}

context=args[1]
bin=paste0(args[2],"bp")
cov=paste0(args[3],"cov")
type=args[4]

print(c(context, bin, cov, type))

library(tidyverse)
library(gplots)

if(bin == '0bp' & cov == '0cov' ){
	files <- dir(pattern=paste(context)) %>% 
		subset(subset=substr(., start=nchar(.)-2, stop=nchar(.)) == 'bed')
} else {
	files <- dir(pattern=paste(context,bin,cov,sep="_"))
}

print(files)

if(type==1){

data <- data_frame(files) %>%
	mutate(file_contents = map(files, read_delim, delim='\t', col_names=F, skip=1)) %>%
	unnest() %>%
	filter(X1 != 'Mt' & X1 != 'ChrM' & X1 != 'Pt' & X1 != 'ChrC') %>%
	select(files, X1, X2, X3, X4) %>%
	mutate(sample=sapply(strsplit(files, '_'), function(l) l[1])) %>%
	mutate(genotype=sapply(strsplit(sample, '-'), function(l) l[1])) %>%
	mutate(rep=sapply(strsplit(sample, '-'), function(l) l[2])) %>%
	mutate(X1=ifelse(substr(X1, start=1, stop=3)=="Chr",paste0(X1),paste0("Chr",X1))) %>%
	na.omit() %>%
	group_by(X1, X2, X3, genotype, rep) %>%
	summarise(met = mean(X4)) %>%
	unite(temp, genotype, rep) %>%
	spread(key=temp, value=met) %>%
	na.omit() %>%
	ungroup() %>%
	select(-X1, -X2, -X3) %>%
	cor() %>%
	as.matrix()

### heatmap
pdf(file=paste0('wig_cor_',context,'_reps.pdf'), width = 0, height = 0, paper="a4r")

heatmap.2(data, 
	trace='none',
	density.info='none',
	symm=F,
	symkey=F,
	key=T,
	colsep = 1:ncol(data),
	rowsep = 1:nrow(data),
	sepcolor = "white",
	sepwidth = c(0.001,0.001),
	dendrogram='both',
	margins = c(8,8),
	cexCol = 1,
	cexRow = 1)
dev.off()

} else {

data <- data_frame(files) %>%
	mutate(file_contents = map(files, read_delim, delim='\t', col_names=F, skip=1)) %>%
	unnest() %>%
	filter(X1 != 'Mt' & X1 != 'ChrM' & X1 != 'Pt' & X1 != 'ChrC') %>%
	select(files, X1, X2, X3, X4) %>%
	mutate(sample=sapply(strsplit(files, '_'), function(l) l[1])) %>%
	mutate(genotype=sapply(strsplit(sample, '-'), function(l) l[1])) %>%
	mutate(rep=sapply(strsplit(sample, '-'), function(l) l[2])) %>%
	mutate(X1=ifelse(substr(X1, start=1, stop=3)=="Chr",paste0(X1),paste0("Chr",X1))) %>%
	na.omit() %>%
	group_by(X1, X2, X3, genotype) %>%
	summarise(met = mean(X4)) %>%
	spread(key=genotype, value=met) %>%
	na.omit() %>%
	ungroup() %>%
	select(-X1, -X2, -X3) %>%
	cor() %>%
	as.matrix()

## heatmap
pdf(file=paste0('wig_cor_',context,'.pdf'), width = 0, height = 0, paper="a4r")

heatmap.2(data,
        trace='none',
        density.info='none',
        symm=F,
        symkey=F,
        key=T,
        colsep = 1:ncol(data),
        rowsep = 1:nrow(data),
        sepcolor = "white",
        sepwidth = c(0.001,0.001),
        dendrogram='both',
        margins = c(8,8),
        cexCol = 1,
        cexRow = 1)
dev.off()
}
