## Scripts for analysing RNA-sequencing

### kallisto_pipe
Pipeline for using Kallisto (Bray et al 2016) for alignment-free transcript quantification. See individual version scripts for running updates and changes.

### subread_pipe
Pipeline for using Subread (more specifically Subjunc) for read alignment for subsequent quantification with featureCounts. See individual script version for running changes and modifications.

### featureCounts
Use featureCounts for gene expression quantification. featureCounts_exon performs the same but at the exon-level rather than gene-level.

### featureCounts_to_edgeR
Template script for differential gene expression testing with edgeR (Robinson et al 2010).

### RNAseq_bam_to_100bpwigs
Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs across annotations of interest.

### rel_expression_plots
Get raw read coverage across features of interest from using WIG files as input (see RNAseq_bam_to_100bp wigs).

### RNAseq_bam_to_bedgraph
Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format, subsequently converting to bigWig files (IGV browsing).

### SUPPA_pipe
SUPPA2 pipeline using event-based analysis to detect alternate splicing and isoform usage.

### split_file
R script required for splitting files appropriately for SUPPA_pipe.

### stringtie_pipe
StingTie2 pipeline that performs reference-guided de novo transcript assembly and quantification.

### stringtie_extract_tpm
Supplementary script for stringtie_pipe to extract TPM based quantification from StringTie2 output.


