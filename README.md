# Molecular Biology Codes 

...from a volcano plot enthusiast.
Welcome to my repository, I hope my codes won't make you inaudibly say 'what ?'. :)

There are several folders in this repository that belong to various attempts at molecular biology analysis.

## RNA-Seq-Workflow
### rnaseq.sh

This folder introduces itself with a rnaseq_pe.sh bash script. This script takes fastq files (paired-end reads) and through a pseudo-pipe, the reads will be aligned to a human genome via HISAT2 alignment tool (1). Then, through SAMTOOLS the sam files will be sorted and converted to bam files (2). This machine-readable alignment files then will be assembled through STRINGTIE and utilizing a python script, count matrices (gene and transcript) will be the outputs (3). There will be no intermediate files left such as bam/sam/gtf type of files. From several GBs of data, we are are now left with count matrices that are can be readily given as an input to DESEQ2 Differential Expression analysis tool in R (4).

The only difference of the rnaseq_se.sh file is that it allows the inputs of fastq files that have single-end reads.

### deseq2_workflow

Being mentioned, the count matrices are used as an input for DESEQ2 Diferential Expression analysis tool which is under Bioconductor project in R. Although it has currently *312* lines or *237* sloc (I was today years old when I learned what sloc meant) it is a pretty straightforward code. First, the required libraries and count matrixes are imported and a metadata is constructed. Before moving on to differential expression analysis itself, two approaches in clustering the samples are being done: pearson correlation matrix combined with cluster-heatmap and PCA (principal component analysis) plot. After the differential expression analysis and log shrinkage (apeglm) method is being utilized and significant genes are picked (5). In order to export genes as tabular indie/unknown filetype "excel", the dataframes are arranged. Heatmap of significant genes, volcano plot and top 20 genes in a dot plot fashion are also drawn. 

## Python-RnaSeq-Post-Analysis

This directory contains Heatmap and PCA codes if you want to further modify and/or want to use Python for these purposes instead of R.

## Alpha-Fold-Interaction-Scores

### Interaction Scoring Script

This Python script first utilizes 'ProDy' library to take the AlphaFold-Multimer predicted structures (.pdb files) and divides the two chains (all of the predicted structures are dimers) and returns pdb files to the directory (6). Then with the 'freesasa' library, SASA scores (accessibility of each residue by a molecule that is approximately in water molecule size) for each residue are calculated for each of the pdb file: The main pdb file (dimer structure) and its two chains (monomers) (7). By comparing the SASA scores of each residue of the chains to their main structure, the buried area of the dimers are calculated. The threshold for SASA score difference was chosen after an extensive trial and error. At the end of this, we have interaction interface (interaction residues) for every dimer combination. PAE (predicted aligned error) is one of the outputs of the AlphaFold and they are stored in the json file as in dictionary-like structure. Then, using intrachain PAE value comparison (again: dimer chains vs monomers) a new score is computed. This is also done with plDDT scores (that are stored in the pdb file itself and a similar code from pDockq paper is used to extract those values) (8) I want to emphasize that this only gets intrachain PAE values (and also plDDT) of the interaction interface. These scores are also stored and visualized with histograms.

In short, this script compares the difference between monomers (predicted as a single chain) and dimer chains (predicted as two chains in multimer form) of intrachain PAE & plDDT scores. These scores score can either be computed from all residues, reflecting a protein-wide difference, or from interaction residues. *I* think this introduces an unbiased scoring system that can be used to assess the probability of two proteins interacting or even being dimerized. 

****

#### References

_1- Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015 Apr;12(4):357-60. doi: 10.1038/nmeth.3317. Epub 2015 Mar 9. PMID: 25751142; PMCID: PMC4655817_

_2- Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352_

_3- Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015 Mar;33(3):290-5. doi: 10.1038/nbt.3122. Epub 2015 Feb 18. PMID: 25690850; PMCID: PMC4643835._

_4- Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8._

_5- Anqi Zhu, Joseph G Ibrahim, Michael I Love, Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences, Bioinformatics, Volume 35, Issue 12, June 2019_

_6- Bakan A, Meireles LM, Bahar I., ProDy: Protein Dynamics Inferred from Theory and Experiments. Bioinformatics 2011 27(11):1575-1577._

_7- Mitternacht S. FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Res. 2016 Feb 18;5:189. doi: 10.12688/f1000research.7931.1. PMID: 26973785; PMCID: PMC4776673._

_8- Bryant, P., Pozzati, G. & Elofsson, A. Improved prediction of protein-protein interactions using AlphaFold2. Nat Commun 13, 1265 (2022). https://doi.org/10.1038/s41467-022-28865-w_
