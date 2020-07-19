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

## file location
data_paths <- file.path("Raw_data/")

## organise files
data_files <- dir(data_paths) %>%
  tibble(filename=.) %>% 
  mutate( dir_contents = map(filename, ~ dir(file.path(data_paths, .)))) %>%
  unnest(cols = c(dir_contents)) %>%
  mutate(test = regexpr(pattern = "Amplification-Data", text = dir_contents)) %>%
  mutate(test2 = regexpr(pattern = "sample-sheet", text = dir_contents))

## set plate name to be analysed = folder name
plate_name <- args[1]
outdir <- args[2]
# plate_name <- unique(data_files$filename)[1] ## testing code
# outdir <- "chipPCR_output" ## testing code

if(file.exists(outdir) == F) {dir.create(outdir)}

dir.create(outdir)

### generate sample information
keyfile <- filter(data_files, filename == plate_name) %>%
  filter(test2 > 0) %>% 
  select(-test, -test2) %>%
  mutate(file_contents=map(filename,~ read_delim(file.path(data_paths,filename,dir_contents),delim=' ',col_names=F))) %>%
  unnest(c(file_contents)) %>% clean_names %>% 
  rename(well_name=x1, sample_name=x2, reaction_volume=x3, well_number=x4) %>%
  mutate(filename = paste(sapply(strsplit(filename, "_sample"), function(l) l[1]), well_number, sep='_')) %>% 
  select(-well_number)

## oragnise amplification data and link to sample information
amp_data <- filter(data_files, filename == plate_name) %>%
  filter(test > 0) %>% select(-test, -test2) %>%
  mutate(file_contents=map(filename,~ read_csv(file.path(data_paths, filename, dir_contents), skip = 36, 
                    col_types = list(col_number(), col_number(), col_character(), col_number(), col_number())))) %>%
  unnest(c(file_contents)) %>% clean_names %>% 
  rename(fluorescence=rn) %>% select(-delta_rn) %>% na.omit %>%
  mutate(plate_name = sapply(strsplit(filename, "_Amp"), function(l) l[1])) %>%
  mutate(filename = paste(sapply(strsplit(filename, "_Amp"), function(l) l[1]), well, sep='_')) %>%
  left_join(keyfile, by="filename") %>%
  select(-dir_contents.x, -dir_contents.y, -filename) %>%
  mutate(sample_name = sapply(strsplit(sample_name, "_rep"), function(l) l[1])) %>%
  na.omit

## subset data to specific plate
a1 <- subset(amp_data, plate_name==plate_name)

## assemble reactions per amplicon
my_list <- list()
for(i in unique(a1$target_name)){
  b <- filter(a1, target_name == i) %>%
    mutate(sample_name = paste(well_name, sample_name, i, sep='_')) %>%
    select(sample_name, cycle, fluorescence) %>%
    spread(sample_name, fluorescence) %>%
    as.data.frame
  
  my_list[[i]] <- b
  
}

### chipPCR analysis to derive Ct
cq <- NULL

pdf(file = paste0(file.path(outdir,plate_name),"_plots.pdf"), paper = "a4")
par(mfrow = c(3,2))
for( i in names(my_list)){
  a <- my_list[[i]]
  
  ### CPP functions to pre-process data (e.g. smooth, normalize, remove background, remove outliers)
  res.CPP <- apply(a[, -1], MARGIN = 2, function(l) {CPP(a[, 1], l,
                                                         method="supsmu",     ## Friedmann supersmoother
                                                         trans=T,             ## background slope correction
                                                         method.reg="lmrob",   ## robust linear reg
                                                         bg.range = c(1,22),   ## set cycle range for background
                                                         method.norm = "none", ## no normalization
                                                         bg.outliers = T      ### remove outliers in background
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
    
    ## test Ct filters
    test <- ifelse(out$max < out$threshold, yes = "bad", no=
                     ifelse( out$FDM < 15 | out$SDM < 15, yes = "bad" , no =
                               ifelse( abs(out$SDC-out$FDM)>2, yes = "bad", no = "Good!")))
    
    cq <- rbind(out, cq)
    col <- rainbow(4)
    
    ## plot
    plot(res, xlab = "Cycle", ylab="Fluorescence", ylim=c(0,ymax), main=ifelse(test == "bad", 
                                                    yes = paste(colnames(res.CPP)[h], "removed!", sep='='), 
                                                    no = colnames(res.CPP)[h]))
    abline(v=summ[2], col=col[2])
    abline(h=threshold, col = 'darkred')
    text(summ["SDM"]-2, ymax/2, paste0("Ct=", round(summ["SDM"], 2)), cex=0.8, col = col[2])
    
    # abline(v=summ, col = col)
    # text(summ["FDM"]-2, max(b$res_cpp_h)/1.2, paste0("FDM~", round(summ["FDM"], 2)), cex=0.8, col = col[1])
    # text(summ["SDM"]-2, max(b$res_cpp_h)/2, paste0("SDM~", round(summ["SDM"], 2)), cex=0.8, col = col[2])
    # text(summ["SDm"]-2, max(b$res_cpp_h)/3, paste0("SDm~", round(summ["SDm"], 2)), cex=0.8, col = col[3])
    # text(summ["SDC"]-2, max(b$res_cpp_h)/4, paste0("SDC~", round(summ["SDC"], 2)), cex=0.8, col = col[4])
    
  }
  
}

dev.off()
par(mfrow=c(1,1))
#######################
## assemble Cts and apply QC filters
input <- as.data.frame(cq) %>% rownames_to_column(var = 'sample') %>% 
  mutate(Ct = ifelse(FDM<=15 | SDM<=15, yes = 40, no = SDM)) %>%
  mutate(Ct = ifelse(abs(SDC-FDM)<=2, yes = Ct, no = 40)) %>%
  mutate(Ct = ifelse(max < threshold, yes = 40, no = Ct))

lbls1 <- sapply(input$sample, function(l) { regexpr(l, pattern = "Internal Control") })
lbls1 <- lbls1[lbls1>1]
lbls2 <- sapply(input$sample, function(l) { regexpr(l, pattern = "N gene") })
lbls2 <- lbls2[lbls2>1]
lbls3 <- sapply(input$sample, function(l) { regexpr(l, pattern = "ORF1ab") })
lbls3 <- lbls3[lbls3>1]
lbls <- unlist(list(lbls1,lbls2,lbls3))

input <- input %>% 
  mutate(probe = substr(sample, start = lbls[match(sample, names(lbls))], stop = nchar(sample))) %>%
  mutate(probe = sapply(strsplit(probe, "_"), function(l) l[1])) %>%
  mutate(sample = substr(sample, start = 1, stop= lbls[match(sample, names(lbls))] -2))

## output Ct plot per sample
pdf(file = paste0(file.path(outdir,plate_name), '_Ctvalues.pdf'))
print(
  ggplot(data=input, aes(x=sample, y=Ct, colour=sample)) +
    geom_jitter() + theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7), legend.position = "none") +
    scale_y_continuous(name = "Ct") +
    scale_x_discrete(name = "Sample") +
    facet_wrap(~probe, ncol = 1)
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
