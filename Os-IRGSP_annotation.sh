#!/bin/bash

## Source GFF files from ENSEMBL Genomes [https://plants.ensembl.org/index.html]
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.61.chr.gff3.gz

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

ens <- gffRead('Oryza_sativa.IRGSP-1.0.61.chr.gff3')

# mRNA annotation - protein coding + canonical isoform only
mrna <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	mutate(biotype=getAttributeField(attributes, 'biotype')) %>%
	mutate(tag=getAttributeField(attributes, 'tag')) %>%
	subset(tag == "Ensembl_canonical" & biotype == "protein_coding") %>%
        select(seqname,start,end,Name,feature,strand) %>%
	unique

write.table(mrna,'Oryza_sativa.IRGSP-1.0.61_mRNA_canonical.bed', sep='\t', row.names=F, col.names=F, quote=F)

mrna <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
        mutate(biotype=getAttributeField(attributes, 'biotype')) %>%
        mutate(tag=getAttributeField(attributes, 'tag'))

exon <- subset(ens, ens$feature == 'exon') %>%
        mutate(Name=getAttributeField(attributes, 'Name')) %>%
        mutate(Gene=sapply(strsplit(Name, "-"),function(l) l[1])) %>%
        mutate(Transcript=sapply(strsplit(Name, "-"), function(l) paste(l[1],l[2], sep='-'))) %>%
	mutate(constitutive=getAttributeField(attributes, 'constitutive')) %>%
	mutate(tag = mrna$tag[match(Transcript, mrna$Name)]) %>%
	mutate(biotype = mrna$biotype[match(Transcript, mrna$Name)]) %>%
	subset(biotype == "protein_coding") %>% ## subset for protein-coding
	subset(tag == "Ensembl_canonical" | constitutive == 1) %>% ## subset for exons in canonical isoform or constitituvely expressed
        select(seqname,start,end,Name,feature,strand) %>%
	unique

write.table(exon,'Oryza_sativa.IRGSP-1.0.61_exon-mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

stop_mrna <- subset(ens, feature == "three_prime_UTR") %>%
        mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(tag = mRNA$tag[match(id, mRNA$Name)]) %>%
	mutate(biotype = mRNA$biotype[match(id, mRNA$Name)]) %>%
	subset(biotype == "protein_coding") %>% ## protein-coding
	subset(tag == "Ensembl_canonical") %>% ## canonical isoform
	mutate(test1 = ifelse(strand == "-", end+1, start-1)) %>%
	mutate(test2 = ifelse(strand == "-", end+3, start-3)) %>%
	mutate(start = ifelse(strand == "-", test1, test2)) %>%
	mutate(end = ifelse(strand == "-", test2, test1)) %>%
	mutate(feature = "stop codon") %>%
        select(seqname, start, end, id, feature, strand)

write.table(stop_mrna, "Oryza_sativa.IRGSP-1.0.61_stop.bed", sep='\t', row.names=F, col.names=F, quote=F)

quit()
n


