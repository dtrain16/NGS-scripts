#!/usr/bin/env Rscript
## Template file to perform DEG calling using for RNAseq data with edgeR

## Installation of edgeR
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('edgeR'))

library(tidyverse)
library(edgeR)
library(scatterplot3d)

## Count files from featureCounts
countFiles <- dir(pattern = ".counts")

## Define sample groups with descriptive labels from filenames
sampleGroups <- sapply(strsplit(countFiles, "-"), function(l) l[1])

## The DGElist object
lbls <- sapply(strsplit(countFiles, "_"), function(l) l[1])
dge <- readDGE(countFiles, columns = c(1,7), group = sampleGroups, label=lbls, skip=1)

# Use sample groups to make design matrix
design <- model.matrix(~0 + sampleGroups)
colnames(design) <- unique(sampleGroups)

############## rRNA filter

## rRNA contamination (ENSEMBL annotation has rRNA genes as "ncRNA")
rRNA <- read.delim("~/scripts/At_rRNA_AGIs.txt", head=F)
rRNA <- as.character(rRNA$V1)

## find rRNAs & count
rRNA.tags <- match(rRNA, rownames(dge$counts))
rRNA_counts <- dge$counts[rRNA.tags, ]
rRNA.rates <- (colSums(rRNA_counts)/dge$samples$lib.size)*100

## rRNA filter
dge$counts <- dge$counts[-rRNA.tags, ]
##############

## Remove organelle transcripts (if applicable) and MIR genes
exc.tags <- rownames(dge$counts)[substr(rownames(dge$counts), start=3, stop=3) == "C" | 
				 substr(rownames(dge$counts), start=3, stop=3) == "M" | 
				 substr(rownames(dge$counts), start=3, stop=3) == "R"] 
exc.tags <- match(exc.tags, rownames(dge$counts))
dge$counts <- dge$counts[-exc.tags, ]

## Abundance filter (CPM > 1 in at least 3 samples)
keep <- rowSums(cpm(dge) > 1) > 3
dge <- dge[keep, ]

## Re-calculate lib size based on retained transcripts
dge$samples$lib.size <- colSums(dge$counts)

## TMM Normalization
dge.tmm <- calcNormFactors(dge, method = "TMM")

## Estimate common, trended and tagwise dispersion
dge.tmm.disp <- estimateDisp(dge.tmm, design, verobse=TRUE, robust=TRUE)

# diagnostic plots
pdf("diagnostic_plots.pdf", paper="a4r")
groups <- unique(as.character(sampleGroups))
n.reps <- length(sampleGroups)/length(groups)
plot(rRNA.rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)", xlab = "Sample", col = rep(rainbow(length(groups)), each = n.reps), lwd = 1.5, type = 'b')
mtext(sampleGroups, at = 1:length(sampleGroups), las=2, cex=0.5)
plotBCV(dge.tmm.disp)
mds <- plotMDS(dge.tmm.disp, ndim=3, dim.plot = c(1,2), col = rep(rainbow(length(groups)), each = n.reps))
s3d <- scatterplot3d(x = mds$cmdscale.out[,1:3], main = "3-dimensional MDS", xlab = "dim 1", ylab = "dim 2", zlab = "dim 3", color = rep(rainbow(length(groups)), each = n.reps), type='h', pch=19, lwd=1.5)
text(s3d$xyz.convert(mds$cmdscale.out[,1:3]), labels=sampleGroups, cex=.75, pos=4)
dev.off()

#####################################
#####################################

##########
## exact tests
##########

## single factor exact tests (pairwise comparisons)
et <- exactTest(dge.tmm.disp, pair=as.character(unique(dge.tmm.disp$samples$group)), dispersion="trended")

# Identify significantlly differentially expressed genes
sigdeg <- decideTestsDGE(et, adjust.method = "fdr", p.value = 0.05, lfc = 1)
summary(sigdeg)

## final DEG (from pairwise et)
# Extract top DEGs based on exactTest
# tt <- topTags(et, adjust.method = "fdr", sort.by="logFC", p.value=0.05)
tt <- topTags(et, adjust.method = "fdr", sort.by="none", p.value=0.05, n=dim(et)[1])
tt <- tt$table[abs(tt$table$logFC) >= 1,]

################
## GLM
###############

fit <- glmQLFit(dge.tmm.disp, design, robust=TRUE, dispersion=dge.tmm.disp$trended.dispersion)

###
pdf("glm_QLfit.pdf", paper="a4r")
plotQLDisp(fit)
dev.off()
###

qlf <- glmQLFTest(fit, contrast=NULL) # tests lfc=0
# OR
qlf <- glmTreat(fit, contrast = NULL, lfc = 1, null = "interval") # tests |lfc| >= 1
tt <-  topTags(qlf, adjust.method = "fdr", sort.by="none", p.value=0.05, n=dim(qlf)[1])
tt <- tt$table

##################################### ##################################

## Output table

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

# read in Araport11 gff3
anno <- gffRead("~/Araport11/Araport11_GFF3_genes_transposons.201606.gff")

# Gene annotation
gene <- subset(anno,anno$feature=='gene') %>%
	mutate(locus=getAttributeField(attributes, 'ID')) %>%
	mutate(description=getAttributeField(attributes, 'Note')) %>%
	mutate(type=getAttributeField(attributes, 'locus_type')) %>%
	mutate(primary=getAttributeField(attributes, 'full_name')) %>%
	mutate(alias=getAttributeField(attributes, 'Alias')) %>%
	select('locus','description','type','primary','alias') %>%
	subset(locus %in% rownames(tt)) %>%
	mutate(logFC=tt$logFC[match(locus, rownames(tt))]) %>%
	mutate(logCPM=tt$logCPM[match(locus, rownames(tt))]) %>%
	mutate(PVal=tt$PValue[match(locus, rownames(tt))]) %>%
	mutate(FDR=tt$FDR[match(locus, rownames(tt))]) %>%
	arrange(desc(logFC))

write.table(gene,paste0('RNAseq_edgeR-exacttest_',groups[1],'vs',groups[2],'.txt'), sep='\t', na='', col.names=T, row.names=F, quote=F)

########### Get logCPM or RPKM for all genes (of interest) manually
## Use Araport11/TAIR10 annotation file to obtain gene lengths
gene.lengths <- anno %>%
	subset(feature == "gene") %>%
	mutate(locus=getAttributeField(attributes, 'ID')) %>%
	filter(locus %in% rownames(dge.tmm.disp$counts)) %>%
	mutate(length = end - start) %>%
	select(locus, length) %>%
	na.omit() 

## calculate CPM by group
logcpm <- cpmByGroup(dge.tmm.disp, prior.count=2, log=TRUE, normalized.lib.sizes=TRUE, dispersion=dge.tmm.disp$trended.dispersion)
rpkm_gr <- rpkmByGroup(dge.tmm.disp, gene.length=gene.lengths$length, dispersion=dge.tmm.disp$trended.dispersion)

