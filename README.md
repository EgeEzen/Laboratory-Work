# Molecular Biology Codes 

...from a volcano plot enthusiast.
Welcome to my repository, I hope my codes won't make you inaudibly say 'what ?'. :)

There are several folders in this repository that belong to various attempts at molecular biology analysis.

## RNA-Seq-Workflow
### rnaseq.sh

This folder introduces itself with a rnaseq_pe.sh bash script. This script takes fastq files (paired-end reads) and through a pseudo-pipe, the reads will be aligned to a human genome via HISAT2 alignment tool. Then, through SAMTOOLS the sam files will be sorted and converted to bam files. This machine-readable alignment files then will be assembled through STRINGTIE and utilizing a python script, count matrices (gene and transcript) will be the outputs. There will be no intermediate files left such as bam/sam/gtf type of files. From several GBs of data, we are are now left with count matrices that are can be readily given as an input to DESEQ2 Differential Expression analysis tool in R.

The only difference of the rnaseq_se.sh file is that it allows the inputs of fastq files that have single-end reads.

### deseq2_workflow

Being mentioned, the count matrices are used as an input for DESEQ2 Diferential Expression analysis tool which is under Bioconductor project in R. Although it has currently *312* lines or *237* sloc (I was today years old when I learned what sloc meant) it is a pretty straightforward code. First, the required libraries and count matrixes are imported and a metadata is constructed. Before moving on to differential expression analysis itself, two approaches in clustering the samples are being done: pearson correlation matrix combined with cluster-heatmap and PCA (principal component analysis) plot. After the differential expression analysis and log shrinkage (apeglm) method is being utilized and significant genes are picked. In order to export genes as tabular indie/unknown filetype "excel", the dataframes are arranged. Heatmap of significant genes, volcano plot and top 20 genes in a dot plot fashion are also drawn. 

## Python-RnaSeq-Post-Analysis

This directory contains Heatmap and PCA codes if you want to further modify and/or want to use Python for these purposes instead of R.
