#!/usr/bin/env Rscript
# See BAM_to_EJC.sh to get input files ("10bp.50.bed")
# Determine transripts with differential occupancy of EJC complex, inferred by abundance of 5'P reads 30-25 nt upstream of the 3' end of an exon.

library(tidyverse)
library(piecewiseSEM)

fls <- dir(pattern = "_10bp.5p.bed")

out <- NULL

for(i in fls){
        # Read in file
        input <- read.delim(i, head=F) %>%
        # Remove reads to plastid and mitochondria
        subset(V1 != 'ChrM' & V1 != 'ChrC' & V1 != 'Mt' & V1 != 'Pt') %>%
        mutate(length = V7 - V6) %>%
        # features at least 50 bp in length
        subset(length > 49) %>%
        # calculate position relative to 3' end of feature
        mutate(pos = ifelse(V10 == "+", V2-V7, V6-V3)) %>%
        subset(pos > -30 & pos < -25) %>%	
	group_by(V1, V2, V3, V8) %>% summarise(val=mean(V4)) %>%
	mutate(sample = sapply(strsplit(i, 'Aligned'), function(l) l[1]))
	
out <- rbind(input,out)

}

test <- group_by(out, sample, V8) %>%
	## sum read end depth across EJC window
	summarise(sum_val = sum(val)) %>%
	spread(sample, sum_val)

### NAs = 0 / no coverage
test[is.na(test)] <- 0
## at least 6 read ends in 2 samples
keep <- rowSums(test[2:5] > 5) > 1
y <- test[keep,]

write.table(y, "ejc_rna_dxo1_counts.txt", sep='\t', quote=F, row.names=F)

y1 <-   gather(y, sample, val, -V8) %>%
  mutate(genotype = sapply(strsplit(sample, "_"), function(l) l[1])) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "dxo1")))

## t.test for differential 5'P read depth at EJC region
out1 <- NULL

for(i in unique(y1$V8)){
	tmp1 <- subset(y1, V8 == i)
	fit <- lm(val ~ genotype, data=tmp1)
	tmp2 <- data.frame(summary(fit)$coefficients)
	coef <- rownames(tmp2)[2] ## coefficient for extracted stats
	t <- tmp2$t.value[2] ## t-value for coeff
	p <- tmp2$Pr...t..[2] ## p-value for coeff

tmp_out <- data.frame(id=i, coef, t, p)	

out1 <- rbind(out1, tmp_out)

}

out2 <- mutate(out1, p.adj = p.adjust(p=p, method='fdr')) %>%
	subset(p.adj < 0.05) %>%
	left_join(y, by=c("id"="V8"))

## output
write.table(out2, "ejc_rna_dxo1_ttest.txt", sep='\t', quote=F, row.names=F)


