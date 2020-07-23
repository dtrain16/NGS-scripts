#!/usr/bin/env Rscript

## Derive Ct values from raw fluorescence curves using chipPCR 5-point stencil 
## to interpolate second derivative max (Ct value) per reaction.

options(echo=T)
args=commandArgs(trailingOnly=T)
print(args)

# require libraries
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(janitor))
suppressPackageStartupMessages(require(chipPCR))
suppressPackageStartupMessages(require(xlsx)) ## install.packages('xlsx') for reading excel files
suppressPackageStartupMessages(library(gtools)) ## instal.packages('gtools')

# Input folder with trailing slash expected
# inputFolder <- args[1]
inputFile <- args[1]
# inputFile <- "VazymeNAT/Raw_data/2020-07-17_Vazyme_Plate1_10ul_384well_QuantStudio12K/V0001CO200717.xlsx"
print(inputFile)

## output directory
outdir <- args[2]
# outdir <- "VazymeNAT/test"

if(file.exists(outdir) == F) {dir.create(outdir)}

## Split the path into a list
split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
pathFolders <- split_path(inputFile)
pathFolders <- rev(pathFolders)
print(pathFolders)
plate_name <- last(pathFolders)
plate_name <- sapply(strsplit(plate_name, ".xlsx"),function(l) l[1])
print(plate_name)
data_paths <- dirname(inputFile)
print(data_paths)

## Organise files
data_files <- dir(data_paths) %>%
  tibble(filename=.) %>% 
  mutate( dir_contents = map(filename, ~ file.path(data_paths, .))) %>%
  unnest(cols = c(dir_contents))

print(data_files)

### generate sample information
keyfile <- read.xlsx(file = inputFile, sheetName = "Sample Setup", 
                     startRow = 37, colIndex = c(1:3,7)) %>%  clean_names
keyfile <- filter(keyfile, target_name != "")

## oragnise amplification data and link to sample information
amp_data <- read.xlsx(file = inputFile, sheetName = "Amplification Data", startRow = 37, 
                      colIndex = c(1:4), colClasses = c("numeric","numeric","character","numeric"))

amp_data <- clean_names(amp_data) %>%
  rename(fluorescence=rn) %>%
  na.omit %>% 
  mutate(well_pos = keyfile$well_position[match(well, keyfile$well)]) %>%
  mutate(sample_name = keyfile$sample_name[match(interaction(well, target_name), interaction(keyfile$well, keyfile$target_name))])

## assemble reactions per amplicon
my_list <- list()
for(i in unique(amp_data$target_name)){
  b <- filter(amp_data, target_name == i) %>%
    mutate(sample_name = paste(well_pos, sample_name, i, sep='_')) %>%
    select(sample_name, cycle, fluorescence) %>%
    spread(sample_name, fluorescence) %>%
    as.data.frame
  
  my_list[[i]] <- b
  
}

### chipPCR analysis to derive Ct
cq <- NULL

pdf(file = file.path(outdir, paste0(plate_name,"_plots.pdf")), paper = "a4")
par(mfrow = c(3,2))
for( i in names(my_list)){
  a <- my_list[[i]]
  
  ### CPP functions to pre-process data (e.g. smooth, normalize, remove background, remove outliers)
  res.CPP <- apply(a[, -1], MARGIN = 2, function(l) {CPP(a[, 1], l,
                                                         method="spline",     ## standard cubic spline smooth
                                                         trans=T,             ## background slope correction
                                                         method.reg="lmrob",   ## robust linear reg
                                                         bg.range = c(1,22),   ## set cycle range for background
                                                         method.norm = "none", ## no normalization
                                                         bg.outliers = F  ## do not remove background outliers
  )[["y.norm"]]})
  
  threshold = quantile(t(a[,-1]))[[1]]
  ymax= round(max(t(a[,-1]))*1.05)
  hist(t(a[,-1]), breaks = 250, xlab = "Fluorescence", main = "Fluorescence distribution")
  abline(v=threshold, col = "red")
  mtext(text = "Threshold", col='red')
  
  ## plot smoothed data
  matplot(res.CPP, type = "l", xlab = "Cycle", ylab = "Fluorescence", main = paste("All samples\n",i))
  
  ## interpolate FDM and SDM from smoothed data and hold Cq values
  for( h in 1:ncol(res.CPP)){
    b <- clean_names(data.frame(cycle = a$cycle, res.CPP[, h]))
    res <- inder(b, smooth.method = NULL)
    summ <- summary(res, print = FALSE)
    out <- as.data.frame(t(summ))
    rownames(out) <- colnames(res.CPP)[h]
    out$max <- max(res[,2])
    out$threshold <- threshold
    out$sample <- rownames(out)
    
    ## test Ct filters
    test <- ifelse(out$max < out$threshold, yes = "bad", no=
                     ifelse( out$FDM < 15 | out$SDM < 15, yes = "bad" , no = "Good!"))
    
    cq <- rbind(out, cq)
    col <- rainbow(4)
    
    ## plot
    plot(res, xlab = "Cycle", ylab="Fluorescence", ylim=c(0,ymax), main=ifelse(test == "bad", 
                                                    yes = paste(colnames(res.CPP)[h], "removed!", sep='='), 
                                                    no = colnames(res.CPP)[h]))
    abline(v=summ[2], col=col[2])
    abline(h=threshold, col = 'darkred')
    text(summ["SDM"]-2, ymax/2, paste0("Ct=", round(summ["SDM"], 2)), cex=0.8, col = col[2])
    
  }
  
}

dev.off()
par(mfrow=c(1,1))

#######################
## assemble Cts and apply QC filters
input <- as.data.frame(cq) %>% 
  mutate(Ct = ifelse(FDM<=15 | SDM<=15, yes = 0, no = SDM)) %>%
  mutate(Ct = ifelse(max < threshold, yes = 0, no = Ct)) %>%
  mutate(target = sapply(strsplit(sample, "_"), function(l) l[3])) %>%
  mutate(sample = sapply(strsplit(sample, "_"), function(l) paste(l[1],l[2], sep = '_'))) %>%
  mutate(well = sapply(strsplit(sample, "_"), function(l) l[1])) %>%
  mutate(well = factor(well, levels = mixedsort(unique(well)))) %>%
  arrange(well)

## output Ct plot per sample
pdf(file = paste0(file.path(outdir,plate_name), '_Ctvalues.pdf'))
print(
  ggplot(data=input, aes(x=well, y=Ct, colour=target)) +
    geom_jitter(size=1) + theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) +
    scale_y_continuous(name = "Ct") +
    scale_x_discrete(name = "Sample")
)
dev.off()

## write out Ct values per sample
tibble('File type'='', 'Vazyme raw CT' = '') %>% 
  write_csv(path = paste0(file.path(outdir,plate_name), '_Ctvalues.csv'))

input %>% 
  mutate(result = "") %>% 
  write_csv(path = paste0(file.path(outdir,plate_name), '_Ctvalues.csv'), col_names = T, append = T)

## clean up workspace
rm(list=ls())
