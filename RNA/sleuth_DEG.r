#!/usr/bin/env Rscript
## Template file to perform DEG calling utilizing bootstrapped transcript quantification using kallisto

## http://www.nature.com/doifinder/10.1038/nmeth.4324
## https://pachterlab.github.io/sleuth/download.html

## installation
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

library(sleuth)
library(tidyverse)

## setup metadata file and add path to .h5 output files
metadata <- read_delim(exp_sample_table.txt, sep='\t', header=T) %>%
	mutate(path = file.path("kallisto_quant", sample, '.h5'))
head(metadata)

## add
