#!/usr/bin/env Rscript

## Derive Ct values from raw fluorescence curves using chipPCR 5-point stencil to interpolate second derivative max.

# options(echo=T)
# args=commandArgs(trailingOnly=T)
# print(args)

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
  unnest %>%
  mutate(test = regexpr(pattern = "Amplification-Data", text = dir_contents)) %>%
  mutate(test2 = regexpr(pattern = "sample-sheet", text = dir_contents))

### generate sample information
keyfile <- filter(data_files, test2 > 0) %>% select(-test, -test2) %>%
  mutate(file_contents=map(filename,~ read_delim(file.path(data_paths, filename, dir_contents),delim=' ',col_names=F))) %>%
  unnest %>% clean_names %>% 
  rename(well=x1, sample_name=x2, reaction_volume=x3, well_number=x4, dilution=x5) %>%
  mutate(filename = paste(sapply(strsplit(filename, "_sample"), function(l) l[1]), well_number, sep='_')) %>% 
  select(-well, -well_number)

## oragnise amplification data and link to sample information
amp_data <- filter(data_files, test > 0) %>% select(-test, -test2) %>%
  mutate(file_contents=map(filename,~ read_csv(file.path(data_paths, filename, dir_contents), skip = 36, 
                    col_types = list(col_number(), col_number(), col_character(), col_number(), col_number())))) %>%
  unnest %>% clean_names %>% 
  rename(fluorescence=rn) %>% select(-delta_rn) %>% na.omit %>%
  mutate(plate_name = sapply(strsplit(filename, "_Amp"), function(l) l[1])) %>%
  mutate(filename = paste(sapply(strsplit(filename, "_Amp"), function(l) l[1]), well, sep='_')) %>%
  left_join(keyfile, by="filename") %>%
  select(-dir_contents.x, -dir_contents.y, -filename) %>%
  mutate(sample_name = sapply(strsplit(sample_name, "_rep"), function(l) l[1])) %>%
  na.omit

## get plate name - can turn into for loop
plate_name = unique(amp_data$plate_name)[1]
a <- subset(amp_data, plate_name==plate_name)

## assemble reactions per amplicon
my_list <- list()
for(i in unique(a$target_name)){
  b <- filter(a, target_name == i) %>%
    mutate(sample_name = paste(sample_name, i, well, sep='_')) %>%
    select(sample_name, cycle, fluorescence) %>%
    spread(sample_name, fluorescence) %>%
    as.data.frame
  
  my_list[[i]] <- b
  
}

### chipPCR analysis to derive Ct

cq <- NULL

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
  hist(t(a[,-1]), breaks = 250, xlab = "Fluorescence", main = "Fluorescence distribution")
  abline(v=threshold, col = "darkred")
  
  ## plot smoothed data
  matplot(res.CPP, type = "l", xlab = "Cycle", ylab = "Fluorescence", main = paste("CPP\n",i))
  
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
    
    ## plot all obtained metrics
    plot(res, 
         main=ifelse(test == "bad", yes = paste(colnames(res.CPP)[h], "removed!", sep='--'), no = colnames(res.CPP)[h]))
    abline(v=summ, col = col)
    text(summ["FDM"]-2, max(b$res_cpp_h)/1.2, paste0("FDM~", round(summ["FDM"], 2)), cex=0.8, col = col[1])
    text(summ["SDM"]-2, max(b$res_cpp_h)/2, paste0("SDM~", round(summ["SDM"], 2)), cex=0.8, col = col[2])
    text(summ["SDm"]-2, max(b$res_cpp_h)/3, paste0("SDm~", round(summ["SDm"], 2)), cex=0.8, col = col[3])
    text(summ["SDC"]-2, max(b$res_cpp_h)/4, paste0("SDC~", round(summ["SDC"], 2)), cex=0.8, col = col[4])
    abline(h=threshold, col = 'darkred')
  }
  
}

## assemble Cts and apply QC filters
input <- as.data.frame(cq) %>% rownames_to_column(var = 'sample') %>% 
  mutate(ct = ifelse(FDM<=15 | SDM<=15, yes = 40, no = SDM)) %>%
  mutate(ct = ifelse(abs(SDC-FDM)<=2, yes = ct, no = 40)) %>%
  mutate(ct = ifelse(max < threshold, yes = 40, no = ct))

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

par(mfrow=c(1,1))

## output Ct plot per sample
pdf(file = paste("chipPCR_output/",plate_name, '.pdf', sep=''))
print(
  ggplot(data=input, aes(x=sample, y=ct, colour=sample)) +
    geom_jitter() + theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7), legend.position = "none") +
    scale_y_continuous(name = "Ct") +
    scale_x_discrete(name = "Sample") +
    facet_wrap(~probe, ncol = 1)
)
dev.off()

## write out Ct values per sample
write_csv(input, path = paste("chipPCR_output/",plate_name, '.csv', sep=''))

