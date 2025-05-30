#!/bin/bash

## Source GFF files from ENSEMBL Genomes [https://plants.ensembl.org/info/website/ftp/index.html]
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gff3.gz
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

ens <- gffRead('Arabidopsis_thaliana.TAIR10.59.gff3')

# Chromosome annotation
chr <- subset(ens,ens$feature=='chromosome') %>%
        select(c('seqname','start','end'))

# Gene annotation
gene <- subset(ens, ens$feature == 'gene') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(gene,'Arabidopsis_thaliana.TAIR10.59_gene.bed', sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(mRNA,'Arabidopsis_thaliana.TAIR10.59_mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation - protein coding + canonical isoform only
mRNA_canon <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	mutate(biotype=getAttributeField(attributes, 'biotype')) %>%
	mutate(tag=getAttributeField(attributes, 'tag')) %>%
	subset(tag == "Ensembl_canonical" & biotype == "protein_coding") %>%
        select(seqname,start,end,Name,feature,strand) %>%
	unique

write.table(mRNA_canon,'Arabidopsis_thaliana.TAIR10.59_mRNA_canonical.bed', sep='\t', row.names=F, col.names=F, quote=F)

## non-coding RNA
ncrna <- subset(ens, feature == "lnc_RNA" | feature == "miRNA" | feature == "ncRNA" | feature == "rRNA" | feature == "snoRNA" | feature == "snRNA" | feature == "tRNA" ) %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(ID=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
        mutate(Isoform=sapply(strsplit(ID, "\\."), function(l) l[2])) %>%
	mutate(Isoform=ifelse(is.na(Isoform) == T, yes = 1, no = Isoform)) %>%
        subset(Isoform==1) %>% ## subset for primary transcript isoform
        select(c('seqname','start','end','ID','feature','strand'))

write.table(ncrna,'Arabidopsis_thaliana.TAIR10.59_ncRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

# all exons based on canonical isoform
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
        mutate(biotype=getAttributeField(attributes, 'biotype')) %>%
        mutate(tag=getAttributeField(attributes, 'tag')) 

exon <- subset(ens, ens$feature == 'exon') %>%
        mutate(Name=getAttributeField(attributes, 'Name')) %>%
        mutate(Transcript=sapply(strsplit(Name, "\\."), function(l) paste(l[1],l[2], sep='.'))) %>%
	mutate(constitutive=getAttributeField(attributes, 'constitutive')) %>%
        mutate(tag = mRNA_canon$tag[match(Transcript, mRNA_canon$Name)]) %>%
	subset(tag == "Ensembl_canonical" | constitutive == 1) %>% ## canonical isoform & constitutively expressed
        select(seqname,start,end,Name,feature,strand)

write.table(exon,'Arabidopsis_thaliana.TAIR10.59_exon.bed', sep='\t', row.names=F, col.names=F, quote=F)

exon1 <- mutate(exon, test = sapply(strsplit(Name, "\\."), function(l) l[3])) %>%
        subset(test == "exon1") %>%
        select(-test)

write.table(exon1,'Arabidopsis_thaliana.TAIR10.59_exon1.bed', sep='\t', row.names=F, col.names=F, quote=F)

# exon from protein-coding mRNAs based on canonical isoform
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	mutate(biotype=getAttributeField(attributes, 'biotype')) %>%
	mutate(tag=getAttributeField(attributes, 'tag'))

cd_exon <- subset(ens, ens$feature == 'exon') %>%
        mutate(Name=getAttributeField(attributes, 'Name')) %>%
	mutate(Parent=getAttributeField(attributes, 'Parent')) %>%
        mutate(Gene=sapply(strsplit(Name, "\\."),function(l) l[1])) %>%
        mutate(Transcript=sapply(strsplit(Name, "\\."), function(l) paste(l[1],l[2], sep='.'))) %>%
	mutate(constitutive=getAttributeField(attributes, 'constitutive')) %>%
	subset(Transcript %in% mRNA$Name) %>% ## subset of exons from annotated mRNAs
	mutate(tag = mRNA$tag[match(Transcript, mRNA$Name)]) %>%
	mutate(biotype = mRNA$biotype[match(Transcript, mRNA$Name)]) %>%
	subset(biotype == "protein_coding") %>% ## subset for protein-coding
	subset(tag == "Ensembl_canonical" | constitutive == 1) %>% ## subset for exons in canonical isoform or constitituvely expressed
        select(seqname,start,end,Name,feature,strand) %>%
	unique

write.table(cd_exon,'Arabidopsis_thaliana.TAIR10.59_exon-mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

## CDS annotation

cds <- subset(ens, feature == "CDS") %>% 
	mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(Isoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
        subset(Isoform==1) %>%
	group_by(seqname, strand, id, feature) %>%
	summarise(start = min(start), end = max(end) ) %>%
	select(seqname, start, end, id, feature, strand) %>%
	arrange(seqname, start)

write.table(cds, "Arabidopsis_thaliana.TAIR10.59_cds.bed", sep='\t', row.names=F, col.names=F, quote=F)

# UTR annotation
utr <- subset(ens, feature == "five_prime_UTR" | feature == "three_prime_UTR") %>%
	mutate(id = getAttributeField(attributes, 'Parent')) %>%
	mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(Isoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
        subset(Isoform==1) %>%
	select(seqname, start, end, strand, id, feature)

write.table(utr, "Arabidopsis_thaliana.TAIR10.59_UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

## separate UTRs

utr_5 <- subset(ens, feature == "five_prime_UTR") %>%
        mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(IsoformIsoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
	subset(Isoform==1) %>%
        select(seqname, start, end, id, feature, strand)

write.table(utr_5, "Arabidopsis_thaliana.TAIR10.59_5UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

utr_3 <- subset(ens, feature == "three_prime_UTR") %>%
        mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(Isoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
        subset(Isoform==1) %>%
        select(seqname, start, end, id, feature, strand)

write.table(utr_3, "Arabidopsis_thaliana.TAIR10.59_3UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

# start and STOP codon for primary mRNA isoform
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
        mutate(ID=getAttributeField(attributes, 'ID')) %>%
        mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
        mutate(Isoform=sapply(strsplit(Name, "\\."), function(l) l[2])) %>%
        subset(Isoform==1)

stop_mrna <- subset(ens, feature == "three_prime_UTR") %>%
        mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
	mutate(Isoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
        subset(Isoform==1) %>%
	subset(id %in% mRNA$Name) %>%
	mutate(test1 = ifelse(strand == "-", end+1, start-1)) %>%
	mutate(test2 = ifelse(strand == "-", end+3, start-3)) %>%
	mutate(start = ifelse(strand == "-", test1, test2)) %>%
	mutate(end = ifelse(strand == "-", test2, test1)) %>%
	mutate(feature = "stop codon") %>%
        select(seqname, start, end, id, feature, strand)

write.table(stop_mrna, "Arabidopsis_thaliana.TAIR10.59_stop.bed", sep='\t', row.names=F, col.names=F, quote=F)

start_mrna <- subset(ens, feature == "five_prime_UTR") %>%
        mutate(id = getAttributeField(attributes, 'Parent')) %>%
        mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
        mutate(Isoform=sapply(strsplit(id, "\\."), function(l) l[2])) %>%
        subset(Isoform==1) %>%
        subset(id %in% mRNA$Name) %>%
        mutate(test1 = ifelse(strand == "-", start-1, end+1)) %>%
        mutate(test2 = ifelse(strand == "-", start-3, end+3)) %>%
        mutate(start = ifelse(strand == "-", test2, test1)) %>%
        mutate(end = ifelse(strand == "-", test1, test2)) %>%
        mutate(feature = "start codon") %>%
        select(seqname, start, end, id, feature, strand)

write.table(start_mrna, "Arabidopsis_thaliana.TAIR10.59_start.bed", sep='\t', row.names=F, col.names=F, quote=F)

quit()
n

########################
## use bedtools getfasta to obtain UTR sequences > 10 bp
####################
bedtools getfasta -bedOut -s -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -bed Arabidopsis_thaliana.TAIR10.59_UTR.bed > TAIR10_UTR_seq.bed

sortBed -i TAIR10_UTR.bed | groupBy -g 5-6 -c 1,2,3,4 -o first,first,first,first | awk '{ print $3,$4,$5,$1,$2,$6 }' OFS='\t' > TAIR10_UTR.sorted.grouped.bed

awk '{ if ($5 == "five_prime_UTR" && $3-$2 > 10) { print } }' TAIR10_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10_UTR.sorted.grouped.5p.bed

awk '{ if ($5 == "three_prime_UTR" && $3-$2 > 10) { print } }' TAIR10_UTR.sorted.grouped.bed | awk '!a[$4]++' > TAIR10_UTR.sorted.grouped.3p.bed

bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -name -s -bed TAIR10_UTR.sorted.grouped.5p.bed -fo TAIR10_5pUTR.fa

bedtools getfasta -fi $HOME/ref_seqs/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -name -s -bed TAIR10_UTR.sorted.grouped.3p.bed -fo TAIR10_3pUTR.fa


