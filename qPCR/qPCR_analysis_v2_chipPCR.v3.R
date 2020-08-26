#!/usr/bin/env Rscript

## Derive Ct values from raw fluorescence curves using chipPCR 5-point stencil 
## to interpolate second derivative max (Ct value) per reaction.
options(echo=T, stringsAsFactors = FALSE, java.parameters = "- Xmx2g")
args=commandArgs(trailingOnly=T)
print(args)

# require libraries
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(janitor))
suppressPackageStartupMessages(require(chipPCR))
suppressPackageStartupMessages(require(xlsx))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(patchwork))

# Input folder with trailing slash expected
inputFile <- args[1]
# inputFile <- "VazymeNAT/Raw_data/2020-07-17_Vazyme_Plate1_10ul_384well_QuantStudio12K/V0001CO200717.xlsx"

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
test <- read.xlsx(file = inputFile, sheetName = "Sample Setup", header = F)
srtrw <- as.numeric(row.names(test[test$X1 == "Well",]))

keyfile <- read.xlsx(file = inputFile, sheetName = "Sample Setup", startRow = srtrw, 
                     colIndex = c(1:3,7), 
                     colClasses = c("character","character","character","character")) %>%  clean_names
keyfile <- filter(keyfile, target_name != "")

###################
## oragnise amplification data and link to sample information
test <- read.xlsx(file = inputFile, sheetName = "Amplification Data", header = F)
srtrw <- as.numeric(row.names(test[test$X1 == "Well",]))

amp_data <- read.xlsx(file = inputFile, sheetName = "Amplification Data", startRow = srtrw, 
                      colIndex = c(1:4), colClasses = c("numeric","numeric","character","numeric"))

amp_data <- clean_names(amp_data) %>%
  rename(fluorescence=rn) %>%
  na.omit %>% 
  mutate(well_pos = keyfile$well_position[match(well, keyfile$well)]) %>%
  mutate(sample_name = keyfile$sample_name[match(interaction(well, target_name), interaction(keyfile$well, keyfile$target_name))])

######################
## assemble reactions per amplicon as input to chipPCR
my_list <- list()
for(i in unique(amp_data$target_name)){
  b <- filter(amp_data, target_name == i) %>%
    mutate(sample_name = paste(well_pos, sample_name, i, sep='_')) %>%
    select(sample_name, cycle, fluorescence) %>%
    spread(sample_name, fluorescence) %>%
    as.data.frame
  
  my_list[[i]] <- b
  
}

##################################
###### chipPCR analysis to derive Ct
################################
cq <- NULL

for( i in names(my_list)){
  a <- my_list[[i]]
  a <- a[-1, ]
  
  ### CPP functions to pre-process data (e.g. smooth, normalize, remove background, remove outliers)
  res.CPP <- apply(a[, -1], MARGIN = 2, function(l) {CPP(a[, 1], l,
                                              method="spline",     ## standard cubic spline smooth
                                              trans=T,             ## background slope correction
                                              method.reg="lmrob",   ## robust linear reg
                                              bg.range = c(0,20),   ## set cycle range for background
                                              method.norm = "none", ## no normalization
                                              bg.outliers = F  ## do not remove background outliers
  )[["y.norm"]]})
  
  threshold = quantile(t(a[,-1]))[[1]]
  
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
    cq <- rbind(out, cq)
  }
}

#######################
## assemble Cts and apply QC filters
input <- as.data.frame(cq) %>% 
  mutate(Ct = ifelse(FDM<=15 | SDM<=15, yes = 50, no = SDM)) %>%
  mutate(Ct = ifelse(max < threshold, yes = 50, no = Ct)) %>%
  mutate(target = sapply(strsplit(sample, "_"), function(l) l[3])) %>%
  mutate(sample = sapply(strsplit(sample, "_"), function(l) paste(l[1],l[2], sep = '_'))) %>%
  mutate(well = sapply(strsplit(sample, "_"), function(l) l[1])) %>%
  mutate(well = factor(well, levels = mixedsort(unique(well)))) %>%
  arrange(well)

############################
### assemble output pdf
###########################

pdf(file = file.path(outdir, paste0(plate_name,"_plots.pdf")), paper = "a4")

## calculated Ct per sample
a1 <- ggplot(data=input, aes(x=well, y=Ct, colour=target)) +
  geom_jitter(size=1) + theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) +
  scale_y_continuous(name = "Ct", breaks = seq(0,50,by = 10)) +
  coord_cartesian(ylim = c(0,50)) +
  scale_x_discrete(name = "Sample") +
  ggtitle("Ct values") +
  theme(legend.position = "bottom")

## raw amplification data
tmp <- mutate(amp_data, sample_name = sapply(strsplit(sample_name, ' '), function(l) l[1])) %>%
  mutate(sample_name = sapply(strsplit(sample_name, '-'), function(l) l[1]))

a2 <- ggplot(data=tmp) + aes(x=cycle, y=fluorescence, colour=sample_name) + theme_bw() + 
  geom_point(size=0.4, alpha=0.4) + 
  geom_vline(xintercept = c(32,38), linetype = 'longdash', colour='black') +
  theme(legend.position = "none") +
  facet_wrap(~target_name, ncol = 2) + 
  ggtitle(paste0(plate_name, " Raw data"))

a1 / a2

par(mfrow = c(3,2))

## Repeat CPP to make plots
for( i in names(my_list)){
  a <- my_list[[i]]
  a <- a[-1, ]
  
  ### CPP functions to pre-process data (e.g. smooth, normalize, remove background, remove outliers)
  res.CPP <- apply(a[, -1], MARGIN = 2, function(l) {CPP(a[, 1], l,
                                        method="spline",     ## standard cubic spline smooth
                                        trans=T,             ## background slope correction
                                        method.reg="lmrob",   ## robust linear reg
                                        bg.range = c(0,20),   ## set cycle range for background
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
    plot(res, xlab = "Cycle", ylab="Fluorescence", ylim=c(0,ymax), 
         main=ifelse(test == "bad", 
                     yes = paste(colnames(res.CPP)[h], "removed!", sep='='),
                     no = colnames(res.CPP)[h]))
    abline(v=summ[2], col=col[2])
    abline(v=c(32,38), col="black")
    abline(h=threshold, col = 'darkred')
    text(summ["SDM"]-2, ymax/2, paste0("Ct=", round(summ["SDM"], 2)), cex=0.8, col = col[2])
    
  }
  
}

par(mfrow = c(1,1))

dev.off()

## write out Ct values per sample
tibble('File type'='', 'Vazyme raw CT' = '') %>% 
  write_csv(path = paste0(file.path(outdir,plate_name), '_Ctvalues.csv'))

DecideTest <- function(internal, n_gene, ORF) {if(internal > 38) {x="Bad reaction"} else 
if(internal>1&internal<=38 & n_gene>38 & ORF>38) {x = "Negative"} else
  if(internal>1&internal<=38 & n_gene>33&n_gene<=38 & ORF>33&ORF<=38){x="Borderline"} else
    if(internal>1&internal<=38 & n_gene>33 & ORF>1&ORF<=33){x="ORF1ab only"} else
      if(internal>1&internal<=38 & n_gene>1&n_gene<=33 & ORF>33){x="N only"} else
        if(internal>1&internal<=38 & n_gene>1&n_gene<=33 & ORF>1&ORF<=33){x="Positive"} 
  else{"Undefined"}
}

output <- select(input, sample, Ct, target, well) %>%
  spread(key=target, value = Ct) %>%
  clean_names() %>%
  mutate(sample = sapply(strsplit(sample,"_"), function(l) l[2])) %>%
  mutate(well = factor(well, levels = mixedsort(unique(well)))) %>%
  arrange(well)

result <- as.matrix(mapply(FUN = DecideTest, internal=output$internal_control, output$n_gene, output$orf1ab))

output <- cbind(output, result)

output %>% 
  write_csv(path = paste0(file.path(outdir,plate_name), '_Ctvalues.csv'), col_names = T, append = T)

## clean up workspace
rm(list=ls())
