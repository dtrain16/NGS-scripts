#!/usr/bin/env Rscript

## Template file for analysing qPCR data
## Takes "complete" output from LinReg to compare N0 (starting transcript abundance) for each amplicon of interest between conditions tested

# require libraries
library(tidyverse)
library(janitor)
library(ggplot2)
library(emmeans)
library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(xlsx)

# specify path to data (.csv) and to key file (attribute 
data_path <- "LinReg_output/"
key <- read.xlsx("Exp542_spreadsheet.xlsx", sheetName = "design keyfile", colIndex = c(1:7)) %>%
	tbl_df()

# set variables
house <- "PP2AA3" # housekeeper(s)
ref_geno <- "WT" # reference genotype
ref_time <- 0 # reference timepoint
N0_type <- "n0_indiv_eff" #n0_indiv_eff or n0_mean_eff
genotype_transcript <- "XRN4" # any transcripts that are being used to genotype (e.g. T-DNA knockouts)
tech_variation_filter <- 0.1 # tech variation filter - to be implemented

# linreg output files
linreg_files <- dir(path = data_path)

# read in data and QC
raw <- data_frame(linreg_files) %>%
  mutate(content = map(linreg_files, ~read_csv(file.path(data_path, .), skip=3))) %>%
  mutate(plate = sapply(strsplit(linreg_files, ".csv"), function(l) l[1])) %>%
  unnest() %>%
  clean_names() %>%
  tbl_df() %>%
  filter(sample_use != "0 0 0" & sample_use != "0 0 3")

# qc raw points (make sure no -999 or 1s & see variation)
ggplot(raw, aes(y=log(n0_indiv_eff), x=amplicon)) +
  geom_boxplot() +
  facet_wrap(~plate) +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  scale_y_continuous(name = "log2 N0")

# summarise tech reps
dat <- mutate(raw, well = sapply(strsplit(as.character(name), " "), function(l) l[1])) %>%
	mutate(name = sapply(strsplit(as.character(name), " "), function(l) l[2])) %>%
	select(plate, well, name, N0_type, amplicon) %>% 
	group_by(plate, amplicon, name) %>%
	summarise(avg = mean(n0_indiv_eff, na.rm=TRUE), 
		  std = sd(n0_indiv_eff, na.rm=TRUE), 
		  n = sum(!is.na(n0_indiv_eff))) %>%
	mutate(genotype = as.factor(key$Genotype[match(name, key$Harvest.ID)])) %>%
	mutate(timepoint = as.factor(key$Time[match(name, key$Harvest.ID)]))

# get reference factor = avg per gene per plate at base comparison (eg WT @ T=0)
ref_factor <- dat %>%
	subset(genotype == ref_geno & timepoint == ref_time) %>%
	group_by(plate, amplicon, genotype, timepoint) %>%
	summarise(ref_factor = mean(avg), 
		  std=sd(avg), 
		  n=sum(!is.na(avg)))

# get normalisation factor based on housekeeper fold-change across conditions
norm_factor <- subset(dat, amplicon == house) %>%
	mutate(ref_factor = ref_factor$ref_factor[match(interaction(plate,amplicon), interaction(ref_factor$plate,ref_factor$amplicon))]) %>%
	mutate(norm = avg/ref_factor) %>%
	group_by(plate, amplicon, genotype, timepoint, ref_factor) %>%
	summarise(norm_factor = mean(norm, na.rm=TRUE), std = sd(norm, na.rm = TRUE), n=sum(!is.na(norm)))

## See tech variation (assuming housekeeper is true)
ggplot(norm_factor, aes(x=timepoint, y=norm_factor, colour=genotype, group=genotype, fill=genotype)) +
	geom_point() + 
	geom_line() + 
	facet_wrap(~plate) + 
	scale_y_continuous(name = paste(house)) + 
	geom_ribbon(aes(ymin=norm_factor - std, ymax=norm_factor + std), alpha=0.2, show.legend = F, linetype='dotted')

## check genotyping transcripts (e.g. 3' downstream of T-DNA)
geno <- filter(dat, amplicon == genotype_transcript) %>%
	mutate(ref_factor = ref_factor$ref_factor[match(interaction(plate,amplicon), interaction(ref_factor$plate, ref_factor$amplicon))]) %>%
	mutate(ref_fc = avg/ref_factor) %>%
	mutate(norm_factor = norm_factor$norm_factor[match(interaction(plate, genotype, timepoint), interaction(norm_factor$plate,norm_factor$genotype,norm_factor$timepoint))]) %>%
	mutate(norm_fc = ref_fc/norm_factor)

ggplot(geno, aes(x=genotype, y=norm_fc, colour=genotype)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~plate + amplicon) + 
  geom_text(aes(label=name), hjust=1.5, size=3)

## calculate normalized fold changes for all target transcripts
to_mdl <- filter(dat, amplicon != house & amplicon != genotype_transcript) %>%
  mutate(ref_factor = ref_factor$ref_factor[match(interaction(plate,amplicon), interaction(ref_factor$plate, ref_factor$amplicon))]) %>%
  mutate(ref_fc = avg/ref_factor) %>%
  mutate(norm_factor = norm_factor$norm_factor[match(interaction(plate, genotype, timepoint), interaction(norm_factor$plate,norm_factor$genotype,norm_factor$timepoint))]) %>%
  mutate(norm_fc = ref_fc/norm_factor) %>%
  mutate(tray = key$Tray[match(name, key$Harvest.ID)]) %>%
  mutate(pos = key$Position[match(name, key$Harvest.ID)]) 
  
ggplot(data = to_mdl, aes(x=timepoint, y=log2(norm_fc), colour=genotype, group=genotype, fill=genotype)) +
      geom_point() +
      facet_wrap(~plate + amplicon, ncol=2) +
      geom_smooth(method = 'loess', formula = y ~ x, linetype = 'dashed', alpha = 0.1)


