#!/usr/bin/env Rscript
# Filter reads from GRAFT-NAD-seq libraries based on presence of A within a defined distance
# Trim filtered reads to the first A at the 5' end of reads

library(tidyverse)
args=commandArgs(trailingOnly=T)
d <- as.numeric(paste(args[1]))
fastq_data <- readLines(args[2])

# Check that the file has a proper number of lines
if (length(fastq_data) %% 4 != 0) {stop("The FASTQ file is not properly formatted.")}

# Group lines into a matrix (4 rows per read)
reads_matrix <- matrix(fastq_data, ncol = 4, byrow = TRUE)

processed_reads <- apply(reads_matrix, 1, function(l) {
	head <- l[1]
	seq <- l[2]
	opt <- l[3]
	qual <- l[4]

	# get sequences that contain 'A" within the first #d bases
	if (grepl("A", substr(seq, 1, d))){
		
		# Find the position of the first 'A'		
		pos <- regexpr("A", seq)[1]
		
		# replace original sequence with trimmed sequence
	        if (pos > 0) {
                
		# Trim sequence to position of first 'A'
                trimmed_sequence <- substr(seq, pos, nchar(seq))
                trimmed_quality <- substr(qual, pos, nchar(qual))
		}
	
		return(c(head,trimmed_sequence,opt,trimmed_quality))} 
	else {return(NULL) }
	}
)

processed_reads <- unlist(processed_reads)

# Write the filtered reads to a new FASTQ file
name <- sapply(strsplit(args[2], ".fastq"), function(l) l[1])
writeLines(processed_reads, paste0(name,"_processed.fq"))
cat(paste0("Filtered and trimmed reads saved to '", name,"_processed.fq'.\n"))

