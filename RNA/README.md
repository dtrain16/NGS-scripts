## Scripts for analysing RNA-sequencing

### kallisto_pipe_v*.sh
Pipeline for using Kallisto (Bray et al 2016) for alignment-free transcript quantification. See individual version scripts for running updates and changes.

### RNAseq_subread_v*.sh
Pipeline for using Subread (more specifically Subjunc) for read alignment for subsequent quantification with featureCounts. See individual script version for running changes and modifications.

### featureCounts.sh
Use featureCounts for gene expression quantification. featureCounts_exon performs the same but at the exon-level rather than gene-level.

### edgeR.r
Template script for differential gene expression testing with edgeR (Robinson et al 2010).

### RNAseq_bam_to_100bpwigs.sh
Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs across annotations of interest.

### rel_expression_plots.r
Get raw read coverage across features of interest from using WIG files as input (see RNAseq_bam_to_100bp wigs).

### RNAseq_bam_to_bedgraph.sh
Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format, subsequently converting to bigWig files (IGV browsing).


