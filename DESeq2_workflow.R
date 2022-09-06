#Installing the DESeq2 and other necessary libraries from Bioconductor project & others

version
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("DESeq2","RColorBrewer", "pheatmap","tidyverse","apeglm","EnhancedVolcano"))
install.packages("ggrepel")

#loading those libraries

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(apeglm)
library(EnhancedVolcano)

#Setting the working directory

setwd("C:/Users/sofiap/R")

#Reading the data
raw_counts <- read.csv("gene_count_matrix_patz1.csv")

#Creating a metadata dataframe
#The row and column names of the dataframe should be changed with respect to experiment

genotype <- c("PATZ1","PATZ1","PATZ1","PATZ1","PATZ1","PATZ1")

condition <- c("WT","WT","WT","KO","KO","KO")

metadata <- data.frame(genotype,condition)

rownames(metadata) <- c("PATZ1_WT1","PATZ1_WT2","PATZ1_WT3",
                               "PATZ1_KO1","PATZ1_KO2","PATZ1_KO3")


#Setting the gene ids of the raw data as the row names (maybe it is different in other datas)
rownames(raw_counts) <- raw_counts$gene_id
raw_counts$gene_id <- NULL

#setting the column names of raw.counts data to row.names of the metadata (this should be changed every analysis
#to fit the experiment) & checking the order of them

colnames(raw_counts) <- c ("PATZ1_KO1" ,"PATZ1_KO2", "PATZ1_KO3","PATZ1_WT1", "PATZ1_WT2", "PATZ1_WT3")
all (colnames(raw_counts) == rownames(metadata))
#this gives FALSE or TRUE

#if FALSE
ind <- match(colnames(raw_counts), rownames(metadata))
metadata <- metadata[ind, ]
all (colnames(raw_counts) == rownames(metadata))
#now it should give TRUE
#be careful about the ORDER of the KOs and WTs. It will be used in further steps.

#if you want to extract these data
#write.csv(as.data.frame(metadata), file = "metadata.csv")
#write.csv(as.data.frame(raw_counts), file = "patz1_rawcounts.csv")

#creating the DESeq2 object
#design is the condition, which is WT vs KO (the base will be controlled in later stages)
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                    colData =  metadata,
                                    design = ~ condition)
#it will give warning message, ignore it (about the base)

#seeing the object
dds

#calculating the median of ratios, and adding to our object 
#(this is already done in while the deseq2 object, but you can use these values separately)
dds <- estimateSizeFactors(dds)

#extracting the normalized counts from object
nrm_counts <- counts(dds, normalized = TRUE)

# Extracting the transformed values (DESEQ2 says blind should be false to reduce noise for downstream analysis)
# DESEQ2 Manual says that bias introduced in blind=FALSE argument is not cruical

vsd <- vst(dds, blind=FALSE) 


#computing pairwise correlation values, first create a matrix from vsd object, then compute correlation
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
#adjust the names if necessary
rownames(vsd_cor) <- c("KO1","KO2","KO3","WT1","WT2","WT3")
colnames(vsd_cor) <- rownames(vsd_cor)

#modifying metadata for only pheatmap use (we will create a hierarchical cluster map from the correlation matrix)
metadata_modified <- metadata
rownames(metadata_modified) <-rownames(vsd_cor)

#then do the correlation heatmap plot
#pheatmap data gets the correlation matrix and in the annotation part it uses the metadata 
#it matches the rownames of the metadataof to column names of the correlation matrix

#select function gets the specific column from a dataframe
pheatmap(vsd_cor, annotation_col = select(metadata_modified, condition),
         main = "Hierarchical Heatmap of PATZ1-KO")
#don't forget the change the name of the map (main argument)

#Doing PCA Analysis
pcaData<-plotPCA(vsd, intgroup="condition",returnData=TRUE)
#when returnData is true, the graph is assigned to a variable, then you can use ggplot to modify it

percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$name <- c("KO1","KO2","KO3","WT1","WT2","WT3")


#now to use ggplot

library(ggrepel)

ggplot(pcaData, aes(PC1, PC2, color=condition,label = name)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(title= "PATZ1 PCA Plot") + geom_label_repel(fontface = "bold", nudge_x = 1 )

#Don't use "color" in the legend, remember to change title of the plot
#ggrepel is very useful if the annotations overlap

#finally doing the DESeq analysis, dds is now the deseq2 object
dds <- DESeq(dds)

#plotting dispersion estimate graph (look further for manual)
plotDispEsts(dds)

#extracting the results
#If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

result <- results(dds, 
                      contrast = c("condition","KO","WT"),
                      alpha = 0.05)
#giving contrast as WT vs KO (WT, the third element should be the base)

#viewing the results
summary(result)
#you can see the up and down regulation here (LFC > 0 is given, we will further limit that)

#Row names are extracted to gene_id column
result$gene_id <- rownames(result)

#extracting this results to an excel file
library(writexl)
#write_xlsx(as.data.frame(result),"patz1_deseq2_genes.xlsx")

#getting the coefficients for the log shrinkage 
resultsNames(dds)
#the second name will be used

#log shrinkage, apeglm method will be used (2018 paper)
result_shr <- lfcShrink(dds,
                         coef="condition_WT_vs_KO",type="apeglm")
#from now on I will continue the analysis with this version!

#investigating the results of the shrinkage version
mcols(result_shr)
summary(result_shr)

#MA plots 
plotMA(result_shr, ylim=c(-2,2))

#plotting counts, if you want to see the change in a specific gene use this function below (change the gene argument)
#plotCounts(dds, gene=which.min(result_shr$svalue), intgroup="condition")

#only getting significants and ordering (FDR 0.05, lfc 1 is picked here)
resOrdered <- result_shr[order(result_shr$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resSig <- subset(resSig, ((log2FoldChange > 1) | (log2FoldChange < -1) ))

summary(resSig)
#88 upregulated
#52 downregulated

#setting the resSig rownames to a gene_names column
resSig$gene_names <- rownames(resSig)
rownames(resSig) <- NULL

#extracting the gene_names column to two: gene_name and gene_id (they were together in my data)
#for example: ENSG00000234711|TUBB8P11 -> ENSG00000234711(in gene_id column) , TUBB8P11 (in gene_name column)

resSig$gene_name <- str_sub(resSig$gene_names,17,-1)
resSig$gene_id <- str_sub(resSig$gene_names, 1, 15)
resSig$gene_names <- NULL

#reordering columns of resSig
colnames(resSig)
resSig <- resSig[, c(7, 6, 1, 2, 3, 4 ,5)]

#seeing the resSig
resSig

#exporting as xlsx
#write_xlsx(as.data.frame(resSig),"patz1_deseq2_genes.xlsx")

#exporting as csv
#write.csv(as.data.frame(resSig), 
#         file="patz1_deseq2_genes.csv")

#getting the normalized&transformed matrix (similar working scheme compared to resSig)
vsd_mat_df <- data.frame(vsd_mat)
vsd_mat_df$gene_names <- rownames(vsd_mat_df)
rownames(vsd_mat_df) <- NULL
vsd_mat_df$gene_name <- str_sub(vsd_mat_df$gene_names,17,-1)
vsd_mat_df$gene_id <- str_sub(vsd_mat_df$gene_names, 1, 15)
vsd_mat_df$gene_names <- NULL
colnames(vsd_mat_df)
vsd_mat_df <- vsd_mat_df[, c(8,7,4,5,6,1,2,3)]

#getting the normalized counts for ONLY significant genes then creating a heatmap

resSig_norm <- merge(as.data.frame(resSig), vsd_mat_df, by ="gene_id")

#setting gene_id to gene_names for the ones that are empty
#for example some genes do not have gene names only gene id: ENSG00000288531, in this case this is set as gene_name
#for loops over rows

for(i in 1:nrow(resSig_norm)) {      
  if (resSig_norm[i ,2 ] ==""){
    resSig_norm[i,2] <- resSig_norm[i,1]
  }
}
#setting gene_name as rownames (going back to there, because I will create volcano plot and heatmap)
rownames(resSig_norm) <- resSig_norm$gene_name.x
rownames(resSig_norm) <- NULL


#we do not need the whole dataset, to use it in a graph, slicing is done

resSig_norm_graph <- resSig_norm[,c(9:14)]

#setting the column names of this new dataframe

colnames(resSig_norm_graph) <- c("KO1","KO2","KO3","WT1","WT2","WT3")

pheatmap(resSig_norm_graph, cluster_rows=TRUE, show_rownames=FALSE,scale = "row",
         cluster_cols=TRUE, annotation_col=select(metadata_modified, condition), 
         main = "PATZ1 Differentially Expressed Genes Heatmap")
#It should be noted that metadata_modified is usede here, similar to heatmap
#show_rownames can be changed to true if you want to annotate them

#getting the lab argument for volcano plot (it requires ALL gene names, same structure as previous two
# maybe a function can be performed)
result_shr_df <- as.data.frame(result_shr) 
result_shr_df$gene_names <- rownames(result_shr_df)
result_shr_df$gene_name <- str_sub(result_shr_df$gene_names,17,-1)
result_shr_df$gene_id <- str_sub(result_shr_df$gene_names, 1, 15)
result_shr_df$gene_names <- NULL

for(i in 1:nrow(result_shr_df)) {      
  if (result_shr_df[i ,7 ] ==""){
    result_shr_df[i,7] <- result_shr_df[i,8]
  }
}


#Doing volcano plot
EnhancedVolcano(result_shr_df,
                lab = result_shr_df$gene_name,
                xlim = c(-10,10),
                pCutoff = 1e-05,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Patz1 Differentially Expressed Genes")

#comparing the results to patz2 (from another set of experiments) 
#VENN diagram will be used
result_shr_df_patz2 <- read.csv("patz2_deg.csv")
library(VennDiagram)

#getting only the gene_ids from both parts and setting to set1 and set2 variables

set1 <- result_shr_df_patz2$gene_id #for patz1
set2 <- resSig_norm$gene_id #for patz2

#we will see the COMMON genes beforehand
lets_look <- merge(result_shr_df_patz2, resSig_norm, by ="gene_id")


#Now to the venn diagram (look further for venn.diagram manual)
#set1 and set2 are used, category.names can be changed and this automatically saves a file (filename argument)
venn.diagram(
  x = list(set1, set2),
  category.names = c("PATZ1 DEGs" , "PATZ2 DEGs"),
  filename="PATZ1vsPATZ2_DEGS.png",
  output=TRUE,
  
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = c("coral1","bisque1"),
  
  # Numbers
  cex = 1.2,
  fontface = "bold",
  fontfamily = "sans",
  
  #set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(27, 27),
  cat.fontfamily = "sans",
)
