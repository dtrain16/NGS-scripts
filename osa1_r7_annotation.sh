#!/bin/bash

## Source GFF files from ENSEMBL Genomes [https://plants.ensembl.org/index.html]
wget https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.gff3.gz

gzip -d *.gz

#### R
library(tidyverse)

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {rv = a[[m]][2]}
         else {rv = as.character(NA)}
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer","integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
        cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

ens <- gffRead('osa1_r7.all_models.gff3')

# mRNA annotation - protein coding + canonical isoform only
mrna <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(Name=getAttributeField(attributes, 'Name')) %>%
        select(seqname,start,end,Name,feature,strand) %>%
	unique

write.table(mrna,'osa1_r7_mRNA_canonical.bed', sep='\t', row.names=F, col.names=F, quote=F)

mrna <- subset(ens, ens$feature == 'mRNA') %>%
	mutate(Name=getAttributeField(attributes, 'Name'))

exon <- subset(ens, ens$feature == 'exon') %>%
	mutate(Name=getAttributeField(attributes, 'ID')) %>%
        mutate(Transcript=sapply(strsplit(Name, ":"),function(l) l[1])) %>%
	mutate(Isoform=sapply(strsplit(Transcript, "\\."),function(l) l[2])) %>%
	subset(Isoform == 1) %>%	
        subset(Transcript %in% mrna$Name) %>% ## exons in mRNA
        select(seqname,start,end,Name,feature,strand) %>%
        unique

write.table(exon,'osa1_r7_exon.bed', sep='\t', row.names=F, col.names=F, quote=F)

stop_mrna <- subset(ens, feature == "three_prime_UTR") %>%
        mutate(Name = getAttributeField(attributes, 'ID')) %>%
        mutate(Transcript = sapply(strsplit(Name, ":"), function(l) l[1])) %>%
        mutate(Isoform=sapply(strsplit(Transcript, "\\."),function(l) l[2])) %>%
        subset(Isoform == 1) %>%
        mutate(test1 = ifelse(strand == "-", end+1, start-1)) %>%
        mutate(test2 = ifelse(strand == "-", end+3, start-3)) %>%
        mutate(start = ifelse(strand == "-", test1, test2)) %>%
        mutate(end = ifelse(strand == "-", test2, test1)) %>%
        mutate(feature = "stop codon") %>%
        select(seqname, start, end, Name, feature, strand)

write.table(stop_mrna, "osa1_r7_stop.bed", sep='\t', row.names=F, col.names=F, quote=F)


quit()
n


