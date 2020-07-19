## convert LC480 output to standardized format for analysis

options(stringsAsFactors = FALSE)
library(tidyverse)
library(janitor)

setwd("C://Users/u4667515/Dropbox/Collab_Projects/Covid19/RawData_RSB_v1/")
outdir <- "../RawData_v2_standardised/"

## probe fluorescence files
fls <- dir(pattern = "Probe_LC480_raw-fluorescence")

for(i in 1:length(fls)){
  a <- read_delim(fls[i], delim = '\t', skip=1, col_names = T) %>% 
    clean_names() %>%
    select(sample_pos, cycle_number, x483_533, x523_568, x558_610) %>%
    rename(well = sample_pos, cycle=cycle_number) %>%
    gather(fluorescence_name, fluorescence, -well, -cycle) %>%
    mutate(fluorescence_name = ifelse(fluorescence_name == "x483_533", yes = "FAM_483-533",
                      no = ifelse(fluorescence_name == "x523_568", yes = "HEX_523-568", no ="Red_558-610")))
  
  write_delim(x = a, path = paste(outdir, fls[i], sep = '/'), delim = '\t')
}

## SYBR fluorescence files
fls <- dir(pattern = "SYBR_LC480_raw-fluorescence")

for(i in 1:length(fls)){
  a <- read_delim(fls[i], delim = '\t', skip=1, col_names = T) %>% 
    clean_names() %>%
    filter(!(cycle_number == 1 & temp != 59.90)) %>% ## remove melt curve values QC = plot(cycle_number ~ temp, a)
    select(sample_pos, cycle_number, x483_533) %>%
    rename(well = sample_pos, cycle=cycle_number) %>%
    gather(fluorescence_name, fluorescence, -well, -cycle) %>%
    mutate(fluorescence_name = ifelse(fluorescence_name == "x483_533", yes = "SYBR_483-533", no=""))

    write_delim(x = a, path = paste(outdir, fls[i], sep = '/'), delim = '\t')
}

## SYBR fluorescence files
fls <- dir(pattern = "sample-sheet")

for(i in 1:length(fls)){
  a <- read_delim(fls[i], delim = '\t', col_names = T) %>% 
    clean_names() %>%
    select(general_pos, general_sample_name, general_filt_comb, general_target_name) %>%
    mutate(general_target_name = ifelse(is.na(general_target_name)==F, yes = general_target_name,
             no= ifelse(general_filt_comb == "483_533", yes = "FAM_RdRP",
                  no = ifelse(general_filt_comb == "523_568", yes = "VIC_human_RP", 
                          no ="ROX_N")))) %>%
    rename(well = general_pos, 
           sample_name=general_sample_name, 
           filter=general_filt_comb, 
           probe_amplicon=general_target_name)
  
  write_delim(x = a, path = paste(outdir, fls[i], sep = '/'), delim = '\t')
}

## view output files to transfer to LabArchives
dir(outdir)

dev.off()
rm(list=ls())

