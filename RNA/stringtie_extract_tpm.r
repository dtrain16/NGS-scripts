## extract TPM values for transcripts assembled by stringtie
## input is stringtie GTF output file
library(tidyverse)
args=commandArgs(trailingOnly=T)

for(i in args){

a <- read.delim(i, header=F, skip=2) %>%
	subset(V3 == "transcript") %>%
	mutate(target_id = sapply(strsplit(as.character(V9), ';'), function(l) l[2])) %>%
	mutate(target_id = sapply(strsplit(target_id, "transcript_id "), function(l) l[2])) %>%
	mutate(length = V5-V4) %>%
	mutate(tpm = sapply(strsplit(as.character(V9), 'TPM '), function(l) l[2])) %>%
	mutate(tpm = as.numeric(sapply(strsplit(tpm, ';'), function(l) l[1]))) %>%
	select(target_id, length, tpm)

outfile <- sapply(strsplit(i, "_stringtie"), function(l) l[1])
write.table(a, paste0(outfile,"_stringtie.tpm"), quote=F, col.names=T, row.names=F)
}

