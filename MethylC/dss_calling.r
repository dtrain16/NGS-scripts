# Script to perform DSS DMR calling 
# Need to enter manually; not setup for running
# Make sure files are converted into right format using DSS_file_prep.r
options(echo=T)
args=commandArgs(trailingOnly=T)
print(args)

# install DSS
# source("http://bioconductor.org/biocLite.R")
# biocLite("DSS")

# Define arguments
context = args[1]
pvalue = args[2]
dlt = args[3]
condition1 = args[4]
condition2 = args[5]

# Read in correctly formatted files
files <- dir(pattern = paste0(context,".output"))

# Define sample groups
group1 <- files[condition1]
group2 <- files[condition2]

library(DSS)
# read input files in DSS format (chr, pos, N, X)
dat1.1 <- read.delim(unlist(group1)[1])
dat1.2 <- read.delim(unlist(group1)[2])
dat1.3 <- read.delim(unlist(group1)[3])

dat2.1 <- read.delim(unlist(group2)[1])
dat2.2 <- read.delim(unlist(group2)[2])
dat2.3 <- read.delim(unlist(group2)[3])

# setup bsseq object
BSobj <- makeBSseqData(list(dat1.1,dat1.2,dat1.3,dat2.1,dat2.2,dat2.3),sampleNames=c("C1","C2","C3","N1","N2","N3"))

# Estimation of methylation means with smoothing by moving averages and smaller smoothing window
dmlTest <- DMLtest(BSobj,group1=c("C1","C2","C3"), group2=c("N1","N2","N3"),smoothing=TRUE,smoothing.span=100)

# identify DMRs based on dmltesting and write out to file
dmrs <- callDMR(dmlTest, delta=dlt, minlen=50, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=pvalue)

## look at distributions of test statistics and p-values
par(mfrow=c(2,2))
hist(dmlTest$stat, 100, main="test statistics")
hist(dmlTest$pval, 100, main="P values")
hist(dmlTest$fdr, 100, main="FDR values")
hist(dmlTest$diff, 100, main="estimates")
dev.off()

# filename
file1=paste0(group1,"vs",group2,"_",context,"_delta=",dlt,"_p=",pvalue,".bed")

# write out file
write.table(dmrs,file=file1,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
