# build annotations based off gff3 files from araport11 release
# https://www.araport.org/data/araport11
# derived from SRE gene_to_gene.sh
# Run lines manually in annotation directory

# Readme file
wget https://www.araport.org/download_file/Araport11_Release_201606/annotation/README.201606.md

# Araport11 annotation in GFF3

# curl -sO -H 'Authorization: Bearer 745dd29759980b058db8fb9efc7af5' https://api.araport.org/files/v2/media/system/araport-public-files//Araport11_Release_201606/annotation/Araport11_GFF3_genes_transposons.201606.gff.gz

wget http://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

gzip -d *gff.gz

# Make bed files

R

library(tidyverse)

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
	cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

ara=gffRead("Araport11_GFF3_genes_transposons.current.gff")

### Gene annotation
gene=subset(ara,ara$feature=='gene')
gene$Name=getAttributeField(gene$attributes, 'Name')
gene$ID=getAttributeField(gene$attributes, 'ID')
gene.out=gene[,c('seqname','start','end','Name','score','strand')]

write.table(gene.out,'Araport11_genes.bed',sep='\t',row.names=F,col.names=F,quote=F)

### TE annotation
te=subset(ara,ara$feature=='transposable_element')
te$Name=getAttributeField(te$attributes, 'Name')
te$ID=getAttributeField(te$attributes, 'ID')
te$alias=getAttributeField(te$attributes, 'Alias')
te$seqname=gsub(pattern="Chr",replacement='', x=te$seqname)
te.out=te[,c('seqname','start','end','Name','score','strand','alias')]

write.table(te.out,'Araport11_TE.bed',sep='\t',row.names=F,col.names=F,quote=F)

### Transcript annotation
mRNA <- subset(ara, feature == "mRNA")
mRNA$name=getAttributeField(mRNA$attributes, 'Name')
mRNA$parent=getAttributeField(mRNA$attributes, 'Parent')
mRNA.out=mRNA[,c('seqname','start','end','name','score','strand')]

write.table(mRNA.out,'Araport11_mRNA.bed',sep='\t',row.names=F,col.names=F,quote=F)

quit()
n

### 5' and 3' UTR annotation
utr <- subset(ara, feature == "five_prime_UTR" | feature == "three_prime_UTR") %>%
	mutate(id = getAttributeField(attributes, 'Parent')) %>%
	select(seqname, start, end, strand, id, feature) %>%
	mutate(id = sapply(strsplit(id, "\\."), function(l) l[1]))
	
	write.table(utr, "Araport11_UTR.bed", sep='\t', row.names=F, col.names=F, quote=F)

## use bedtools getfasta to obtain sequences in utr intervals 
## bedtools getfasta -bedOut -s -fi TAIR10_Chr.all.fasta -bed Araport11_UTR.sorted.bed > Araport11_UTR_seq.bed

##########

rm *gff

