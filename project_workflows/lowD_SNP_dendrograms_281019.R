## check out https://github.com/borevitzlab/brachy-genotyping-notes/blob/master/snprelate.Rmd

library(SNPRelate)
library(tidyverse)

metadata = read.csv("brachy-metadata.csv")
lowd_names = read.delim("lowD_GBS_sample_names.txt")

#snpgdsVCF2GDS("freebayes~GBS~lowD.sorted.vcf.gz",
#              "freebayes~GBS~lowD.gds")

geno = snpgdsOpen("freebayes~GBS~lowD.gds", allow.duplicate = T, readonly = T)
samp = snpgdsSummary(geno)$sample.id
## Filter out very bad missing data
snps = snpgdsSelectSNP(geno, missing.rate=0.999, autosome.only=F)

## Functions for further filtering
ssp.filt = function(geno, samps, snps,  max.snp.miss.rate=0.99,
                    max.samp.miss.rate=0.99, min.maf=0.001 ) {
    miss.samp = snpgdsSampMissRate(geno, snp.id=snps, sample.id=samps)
    
    hist(miss.samp, breaks=100, main="Sample Missing Data (pre-filt)")
    abline(v=max.samp.miss.rate, col="blue", lwd=2)
    
    samps = samps[miss.samp <= max.samp.miss.rate]
    
    srf = snpgdsSNPRateFreq(geno, sample.id=samps, snp.id = snps)
    miss.snp = srf$MissingRate
    hist(miss.snp, breaks=100, main="SNP missing data")
    abline(v=max.snp.miss.rate, col="blue", lwd=2)
    
    maf = srf$MinorFreq
    hist(maf, breaks=50, main="SNP MAF")
    abline(v=min.maf, lwd=2, col="blue")
    
    snps = snpgdsSelectSNP(geno, sample.id=samps, snp.id=snps, maf=min.maf, missing.rate=max.snp.miss.rate, autosome.only=F)
    
    miss.samp = snpgdsSampMissRate(geno, snp.id=snps, sample.id=samps)
    hist(miss.samp, breaks=100, main="Sample Missing Data (post-filt)")
    
    print(paste("Num SNPs:", length(snps)))
    print(paste("Num Samples:", length(samps)))
    return(list(snps=snps, samps=samps, miss.samp=miss.samp))
}

ssp.geno = function(geno, filt) {
    ibs = snpgdsIBS(geno, sample.id=filt$samps, snp.id=filt$snps, autosome.only=F, num.thread=4)
    ibs.nacnt = rowSums(is.na(ibs$ibs))
    table(ibs.nacnt)
    return(ibs)
}

filt.dis = ssp.filt(geno, samp, snps, min.maf=0.01, max.samp.miss.rate = 0.995, max.snp.miss.rate=0.95)
dev.off()

dist <- snpgdsDiss(geno, sample.id=filt.dis$samps, snp.id=filt.dis$snps, autosome.only=F)
dist$sample.id <- paste(lowd_names$acc[match(dist$sample.id, lowd_names$anon)])
hc.dis.plt = snpgdsHCluster(dist) %>% snpgdsCutTree(label.H=F, label.Z=F)

pdf("Brachy_SNP_dendro.pdf", pointsize=4)
snpgdsDrawTree(hc.dis.plt, leaflab="perpendicular", cex.lab=0.01)
dev.off()

