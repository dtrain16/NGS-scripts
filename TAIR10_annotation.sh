#!/bin/bash

## Source GFF files from ENSEMBL Genomes [https://plants.ensembl.org/info/website/ftp/index.html]
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

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

ens <- gffRead('Arabidopsis_thaliana.TAIR10.54.gff3')

# Chromosome annotation
chr <- subset(ens,ens$feature=='chromosome') %>%
        select(c('seqname','start','end'))

# Gene annotation
gene <- subset(ens, ens$feature == 'gene') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(gene,'Arabidopsis_thaliana.TAIR10.54_gene.bed', sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(mRNA,'Arabidopsis_thaliana.TAIR10.54_mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

# exon annotation based on primary isoform
exon <- subset(ens, ens$feature == 'exon') %>%
        mutate(Name=getAttributeField(attributes, 'Name')) %>%
        mutate(Parent=getAttributeField(attributes, 'Parent')) %>%
	mutate(Isoform=sapply(strsplit(Parent, "\\."), function(l) l[2])) %>%
	subset(Isoform==1) %>% ## subset for primary transcript isoform
        select(c('seqname','start','end','Name','score','strand'))

write.table(exon,'Arabidopsis_thaliana.TAIR10.54_exon.bed', sep='\t', row.names=F, col.names=F, quote=F)

exon1 <- mutate(exon, test = sapply(strsplit(Name, "\\."), function(l) l[3])) %>%
	subset(test == "exon1") %>%
	select(-test)

write.table(exon1,'Arabidopsis_thaliana.TAIR10.54_exon1.bed', sep='\t', row.names=F, col.names=F, quote=F)

# UTR annotation
utr <- subset(ens, feature == "five_prime_UTR" | feature == "three_prime_UTR") %>%
	mutate(id = getAttributeField(attributes, 'Parent')) %>%
	mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	select(seqname, start, end, strand, id, feature)

write.table(utr, "Arabidopsis_thaliana.TAIR10.54_UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

quit()
n

########################
## use bedtools getfasta to obtain UTR sequences > 10 bp
####################
bedtools getfasta -bedOut -s -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -bed Arabidopsis_thaliana.TAIR10.54_UTR.bed > TAIR10_UTR_seq.bed

sortBed -i TAIR10_UTR.bed | groupBy -g 5-6 -c 1,2,3,4 -o first,first,first,first | awk '{ print $3,$4,$5,$1,$2,$6 }' OFS='\t' > TAIR10_UTR.sorted.grouped.bed

awk '{ if ($5 == "five_prime_UTR" && $3-$2 > 10) { print } }' TAIR10_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10_UTR.sorted.grouped.5p.bed

awk '{ if ($5 == "three_prime_UTR" && $3-$2 > 10) { print } }' TAIR10_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10_UTR.sorted.grouped.3p.bed

bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -name -s -bed TAIR10_UTR.sorted.grouped.5p.bed -fo TAIR10_5pUTR.fa

bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -name -s -bed TAIR10_UTR.sorted.grouped.3p.bed -fo TAIR10_3pUTR.fa

# clean up
rm TAIR10_UTR.sorted.*.bed -v

