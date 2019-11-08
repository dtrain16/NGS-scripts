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
    
    snps = snpgdsSelectSNP(geno, sample.id=samps, snp.id=snps, maf=min.maf,
                                     missing.rate=max.snp.miss.rate, autosome.only=F)
    
    miss.samp = snpgdsSampMissRate(geno, snp.id=snps, sample.id=samps)
    hist(miss.samp, breaks=100, main="Sample Missing Data (post-filt)")
    
    print(paste("Num SNPs:", length(snps)))
    print(paste("Num Samples:", length(samps)))
    return(list(snps=snps, samps=samps, miss.samp=miss.samp))
}

ssp.geno = function(geno, filt) {
    ibs = snpgdsIBS(geno, sample.id=filt$samps, snp.id=filt$snps,
                    autosome.only=F, num.thread=4)
    ibs.nacnt = rowSums(is.na(ibs$ibs))
    table(ibs.nacnt)
    return(ibs)
}

distimpute = function(ibs, thresh=5, maxdist = 1) {
    dist = 1 - ibs$ibs 
    nasum = colSums(is.na(dist))
    pass = nasum < thresh
    dist = dist[pass,pass]
    
    ibs$sample.id = ibs$sample.id[pass]
    
    dist.ut = dist
    dist.ut[upper.tri(dist.ut,diag=T)] = 0
    nasum = colSums(is.na(dist.ut))
    
    # impute
    for (j in which(nasum > 0)) {
        for (i in which(is.na(dist[,j]))) {
            k = which(dist.ut[,j] == max(dist.ut[,j], na.rm=T))
            if (length(k) > 1) {
                k = k[1]
            }
            if (dist[k, j] <= maxdist) {
                dist[i, j] = max(dist[k, j], # k, j is neighbour -> self
                                 dist[k, i]) # k, i is neighbour -> other
                dist[j, i] = max(dist[k, j], # k, j is neighbour -> self
                                 dist[k, i]) # k, i is neighbour -> other
            }
        }
    }
    
    ibs$ibs = 1 - dist
    return(ibs)
}


distimpute2 = function(ibs, max.NAs=0, max.dist = 0.2) {
    dist = 1 - ibs$ibs 
    
    dist.ut = dist
    dist.ut[upper.tri(dist.ut,diag=T)] = 0
    nasum = colSums(is.na(dist.ut))
    num.imputed = 0
    # impute
    for (j in which(nasum > 0)) {
        for (i in which(is.na(dist[,j]))) {
            k = which(dist.ut[,j] == max(dist.ut[,j], na.rm=T))
            if (length(k) > 1) {
                k = k[1]
            }
            if (dist[k, j] <= max.dist) {
                num.imputed = num.imputed + 1
                dist[i, j] = max(dist[k, j], # k, j is neighbour -> self
                                 dist[k, i]) # k, i is neighbour -> other
                dist[j, i] = max(dist[k, j], # k, j is neighbour -> self
                                 dist[k, i]) # k, i is neighbour -> other
            }
        }
    }

    num.removed = 0
    while (sum(is.na(dist)) > 0) {
        rm = which.max(colSums(is.na(dist)))
        dist = dist[-rm,]
        dist = dist[,-rm]
        ibs$sample.id = ibs$sample.id[-rm]
        num.removed = num.removed + 1
    }
    ibs$ibs = 1 - dist
    ibs$num.imputed = num.imputed
    ibs$num.removed = num.removed
    print(paste("Num imputed:", num.imputed))
    print(paste("Num removed:", num.removed))
    return(ibs)
}

filt.dis = ssp.filt(geno, samp, snps, min.maf=0.01, max.samp.miss.rate = 0.995, 
			max.snp.miss.rate=0.95)

ibs.dis = ssp.geno(geno, filt.dis)
ibs.dis.imp = distimpute2(ibs.dis, max.dist = 0.2)

ibs.dis.plt = ibs.dis.imp
ibs.dis.plt$sample.id <- paste(lowd_names$acc[match(ibs.dis.imp$sample.id, lowd_names$anon)])
hc.dis.plt = snpgdsHCluster(ibs.dis.plt) %>% snpgdsCutTree(z.threshold=20, label.H=F, label.Z=F)

pdf("SNP-dendro.pdf")
snpgdsDrawTree(hc.dis.plt, leaflab="perpendicular")
dev.off()

