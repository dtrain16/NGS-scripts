#!/usr/bin/env Rscript
# Filter reads from GRAFT-NAD-seq libraries based on presence of A within a defined distance
# TO BE IMPLEMENTED: Trim filtered reads to the first A at the 5' end of reads

library(tidyverse)
args=commandArgs(trailingOnly=T)
d <- as.numeric(paste(args[1]))

fastq_data <- readLines(args[2])

# Check that the file has a proper number of lines
if (length(fastq_data) %% 4 != 0) {stop("The FASTQ file is not properly formatted.")}

# Create a container for the filtered reads
filtered_reads <- NULL

# Loop through the reads (FASTQ lines are in groups of 4)
for (i in seq(from = 1, to = length(fastq_data), by = 4)) {
 
	# Extract sequence line (2nd line in the block)
	sequence <- fastq_data[i + 1]
       
	# Check if the first five bases (after adapter) do contain 'A'
	if (grepl("A", substr(sequence, 1, d))) {
        
		# Append all 4 lines (identifier, sequence, '+', and quality) to filtered_reads
        	filtered_reads <- c(filtered_reads, fastq_data[i:(i + 3)])
		}
}

# Create a container for the trimmed reads
trimmed_reads <- NULL

# Loop through trimmed reads (in groups of 4) and trim to A
for (i in seq(from = 1, to = length(filtered_reads), by = 4)) {

	# Extract sequence line (2nd line in the block)
        sequence <- filtered_reads[i + 1]
	
	# Extract quality scores	
	quality <- filtered_reads[i+3]

	# Find the position of the first 'A'
	pos <- regexpr("A", sequence)[1]

	# replace original sequence with trimmed sequence	
	if (pos > 0) {
		# Trim sequence to position of first 'A'
		trimmed_sequence <- substr(sequence, pos, nchar(sequence))
		trimmed_quality <- substr(quality, pos, nchar(quality))		

		#  Append all 4 lines (identifier, sequence, '+', and quality) to trimmed_reads
		trimmed_reads <- c(trimmed_reads, filtered_reads[i], trimmed_sequence, filtered_reads[i+2], trimmed_quality ) 
	}
}

# Write the filtered reads to a new FASTQ file
#writeLines(filtered_reads, paste0(args[2],"_filtered.fq")
#cat(paste0("Filtered reads saved to '", args[2],"_filtered.fq'.\n"))

# Write the filtered reads to a new FASTQ file
name <- sapply(strsplit(args[2], ".fq"), function(l) l[1])
writeLines(trimmed_reads, "trimmed_output.fq")
cat(paste0("Filtered and trimmed reads saved to '", name,"_trimmed.fq'.\n"))

