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

## Alpha-Fold-Interaction-Scores

### Interaction Scoring Script

This Python script first utilizes 'ProDy' library to take the AlphaFold-Multimer predicted structures (.pdb files) and divides the two chains (all of the predicted structures are dimers) and returns pdb files to the directory. Then with the 'freesasa' library, SASA scores (accessibility of each residue by a molecule that is approximately in water molecule size) for each residue are calculated for each of the pdb file: The main pdb file (dimer structure) and its two chains (monomers). By comparing the SASA scores of each residue of the chains to their main structure, the buried area of the dimers are calculated. The threshold for SASA score difference was chosen after an extensive trial and error. PAE (predicted aligned error) is one of the outputs of the AlphaFold and they are stored in the json file as in dictionary-like structure. Making use of these scores, an interaction score between the predicted dimers is formed. It should be noted that for PAE, only the residues that are found as interaction sites are chosen (from buried area,SASA) and the PAE scores were computed from interchain residues. Histograms of said scores are also constructed.

In short, this script finds the interaction residues of dimer predicted structures (AlphaFold-Multimer predictions) and computes an 'interaction score' from interaction interface interchain PAE comparison, aimed to reflect the probability of two proteins to form a dimer. (I know that this sentence is 'inter'crazy).  

### Monomer Dimer Comparison

This script is an evolved version of the previous one, because the results of the previous scores were not truly representative due to better prediction of some structures of AlphaFold. I aimed to be freed of this bias by making a new approach. This time, the monomer structures of the proteins was predicted by AlphaFold singularly and the dimer chains (which are monomers from dimer structures) were compared to their singular forms. Before going into PAE value calculations, RMSD (root mean square deviation) values were computed (dimer chains vs monomers) to reflect how first and second chains are going under structural changes to fit into each other. Then, using intrachain PAE value comparison (again: dimer chains vs monomers) a new score is computed. I want to emphasize that this score does not only utilize interaction residues but all residues of the protein. There is however an algorithm in the later stages of the script that only gets intrachain PAE values of the intearction interface (those residues are extracted from the first script). These scores are also stored and visualized with histograms.

In short, this script compares the difference between monomers (predicted as a single chain) and dimer chains (predicted as two chains in multimer form) of intrachain PAE scores. This PAE score can either be computed from all residues, reflecting a protein-wide difference, or from interaction residues that were established in the previous script. *I* think this introduces an unbiased scoring system that can be used to assess the probability of two proteins interacting or even being dimerized. 
