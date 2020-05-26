#!/usr/bin/env Rscript

## Template R code for calculating Ct values from raw fluorescence curves using chipPCR 5-point stencil to interpolate second derivative max.

# require libraries
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(janitor))
suppressPackageStartupMessages(require(chipPCR))

## set working directory
workdir <- "./"
setwd(workdir)




