#!/bin/bash

## Source GFF files from ENSEMBL Genomes
wget ftp://ftp.ensemblgenomes.org/pub/release-44/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz

gzip -d Arabidopsis_thaliana.TAIR10.44.gff3.gz

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

ens <- gffRead('Arabidopsis_thaliana.TAIR10.44.gff3')

# Chromosome annotation
chr <- subset(ens,ens$feature=='chromosome') %>%
        select(c('seqname','start','end'))

# Gene annotation
gene <- subset(ens, ens$feature == 'gene') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(gene,'TAIR10.44_genes.bed', sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(mRNA,'TAIR10.44_mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

# UTR annotation
utr <- subset(ens, feature == "five_prime_UTR" | feature == "three_prime_UTR") %>%
	mutate(id = getAttributeField(attributes, 'Parent')) %>%
	mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	select(seqname, start, end, strand, id, feature)

write.table(utr, "TAIR10.44_UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

quit()
n

########################
## use bedtools getfasta to obtain UTR sequences > 10 bp
####################
sortBed -i TAIR10.44_UTR.bed | groupBy -g 5-6 -c 1,2,3,4 -o first,first,first,first | awk '{ print $3,$4,$5,$1,$2,$6 }' OFS='\t' > TAIR10.44_UTR.sorted.grouped.bed
awk '{ if ($5 == "five_prime_UTR" && $3-$2 > 10) { print } }' TAIR10.44_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10.44_UTR.sorted.grouped.5p.bed
awk '{ if ($5 == "three_prime_UTR" && $3-$2 > 10) { print } }' TAIR10.44_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10.44_UTR.sorted.grouped.3p.bed
bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Athal.TAIR10.44.dna.fa -name -s -bed TAIR10.44_UTR.sorted.grouped.5p.bed -fo TAIR10.44_5pUTR.fa
bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Athal.TAIR10.44.dna.fa -name -s -bed TAIR10.44_UTR.sorted.grouped.3p.bed -fo TAIR10.44_3pUTR.fa

# clean up
rm *gff3 TAIR10.44_UTR.sorted.*.bed -v

