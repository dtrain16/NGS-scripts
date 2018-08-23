#!/usr/bin/env Rscript
# Produce DSS input files from BED files
# run in folder with cov files of interest and tell which context you want to merge together
# USAGE: DSS_file_prep.r <context>

options(echo=T)
args=commandArgs(trailingOnly=T)
print(args)

## mC context to test
context = args[1]

## get cov files 
files=dir(pattern=paste0(context,".bed"))

## get total and met counts for each sample in given context and write out as separate files
for(i in 1:length(files)){
file <- read.delim(files[i], head=F)
file <- file[file$V1 != "ChrM",]
file <- file[file$V1 != "ChrC",]
file <- file[,c(1,2,7,5)]
test <- as.numeric(regexec(text = as.character(files[i]), pattern=".bed"))
sample <- substr(as.character(files[i]), start = 1, stop = test-1)
colnames(file)=c('chr','pos','N','X')
write.table(x=file, file=paste0(sample,"_output.txt"),sep='\t', quote = F, col.names=T, row.names=F)
}
