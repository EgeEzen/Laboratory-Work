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
library(writexl)
library(ggrepel)

#Setting the working directory

setwd("C:/Users/sofiap/R/sutlulab")

#Reading the data
raw_counts_1 <- read.csv("gene_count_matrix_sutlulab.csv")
colnames(raw_counts_1) <- c ("gene_id","WT_1","DGEZ_1","iT2puro_1","DGE_TYR_1","DGEZ_TYR_1","DGEZ_2","DGE_1","DGEZ_TYR_2","iT2puro_2","DGE_TYR_2","WT_2","WT_3","DGEZ_TYR_3","DGE_2","DGEZ_3","DGE_TYR_3")
raw_counts_1 <- raw_counts_1[,c(1,2,12,13,3,7,16)]

raw_counts_2 <- read.csv("gene_count_matrix_sutlulab_1.csv")
colnames(raw_counts_2) <- c("gene_id" , "iT2puro_3","DGE_3")
raw_counts_2 <- raw_counts_2[,c(1,3)]

raw_counts <- raw_counts_1
rownames(raw_counts) <- raw_counts$gene_id  
raw_counts$gene_id <- NULL

#Creating a metadata dataframe
#The row and column names of the dataframe should be changed with respect to experiment

genotype <- c("NK92","NK92","NK92","NK92","NK92","NK92")
condition <-  c("WT","WT","WT","DGEZ","DGEZ","DGEZ")
metadata <- data.frame(genotype,condition)
rownames(metadata) <- c ("WT_1", "WT_2","WT_3","DGEZ_1","DGEZ_2","DGEZ_3")
#setting the column names of raw.counts data to row.names of the metadata (this should be changed every analysis
#to fit the experiment) & checking the order of them

all (colnames(raw_counts) == rownames(metadata))
#now it should give TRUE
#be careful about the ORDER of the KOs and WTs. It will be used in further steps.

#creating the DESeq2 object
#design is the condition, which is WT vs KO (the base will be controlled in later stage)
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                    colData =  metadata,
                                    design = ~ condition)
#it will give warning message, ignore it (about the base)

#seeing the object
dds

#calculating the median of ratios, and adding to our object 
#(this is already done in while the deseq2 object, but you can use these values separately)
dds <- estimateSizeFactors(dds)

#factor level is set here, this means that the base is now untreated (WT) condition
dds$condition <- relevel(dds$condition, ref = "WT")


#extracting the normalized counts from object
nrm_counts <- counts(dds, normalized = TRUE)

# Extracting the transformed values (DESEQ2 says blind should be false to reduce noise for downstream analysis)
# DESEQ2 Manual says that bias introduced in blind=FALSE argument is not cruical

vsd <- vst(dds, blind=FALSE) 


#computing pairwise correlation values, first create a matrix from vsd object, then compute correlation
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
#adjust the names if necessary


#then do the correlation heatmap plot
#pheatmap data gets the correlation matrix and in the annotation part it uses the metadata 
#it matches the rownames of the metadataof to column names of the correlation matrix

#select function gets the specific column from a dataframe
pheatmap(vsd_cor, annotation_col = select(metadata, condition),
         main = "Hierarchical Heatmap of from Wildtype vs DGEZ", filename ="Wildtype_vs_DGEZ_Samples_clustermap.tiff",width= 10, height = 8)
#don't forget the change the name of the map (main argument)

#Doing PCA Analysis
pcaData<-plotPCA(vsd, intgroup="condition",returnData=TRUE)
#when returnData is true, the graph is assigned to a variable, then you can use ggplot to modify it

percentVar <- round(100 * attr(pcaData, "percentVar"))

#now to use ggplot
jpeg(file="WT_vs_DGEZ_PCA_plot.tiff",width = 6000, height= 6000,res = 600)

ggplot(pcaData, aes(PC1, PC2, color=condition,label = name)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(title= "Wild Type vs DGEZ PCA Plot") +
  geom_label_repel(fontface = "bold", nudge_x = 1, show.legend = FALSE )
dev.off()
#Don't use "color" in the legend, remember to change title of the plot
#ggrepel is very useful if the annotations overlap

#finally doing the DESeq analysis, dds is now the deseq2 object
dds <- DESeq(dds)

#plotting dispersion estimate graph (look further for manual)
plotDispEsts(dds)

#extracting the results
#If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

result <- results(dds, alpha = 0.05)

#viewing the results
summary(result)
result

#you can see the up and down regulation here (LFC > 0 is given, we will further limit that)

#Row names are extracted to gene_id column
result$gene_id <- rownames(result)

#extracting this results to an excel file

#write_xlsx(as.data.frame(result),"patz1_deseq2_genes.xlsx")

#getting the coefficients for the log shrinkage 
resultsNames(dds)

#the second name will be used

#log shrinkage, apeglm method will be used (2018 paper)
result_shr <- lfcShrink(dds,
                         coef="condition_DGEZ_vs_WT",type="apeglm")
#from now on I will continue the analysis with this version!

#investigating the results of the shrinkage version (please check if it is condition treated vs untreated)
result_shr
mcols(result_shr)
summary(result_shr)

#MA plots 
plotMA(result_shr, ylim=c(-2,2))

#plotting counts, if you want to see the change in a specific gene use this function below (change the gene argument)
#plotCounts(dds, gene=which.min(result_shr$svalue), intgroup="condition")

#only getting significants and ordering (FDR 0.05, lfc 1 is picked here)
resSig <- result_shr[order(result_shr$padj),]
resSig <- subset(resSig, padj < 0.05)
resSig <- subset(resSig, ((log2FoldChange > 1) | (log2FoldChange < -1) ))

summary(resSig)
#57 upregulated
#84 downregulated

#a dataframe with genenames as its rownames and non-organized gene_id names like ENSG00000234711|TUBB8P11
data_frame_me <- function(df){
  df_column_number <- length(colnames(df))
  df <- data.frame(df)
  df$gene_names <- rownames(df)#setting the rownames to a gene_names column
  rownames(df) <- NULL
  #extracting the gene_names column to two: gene_name and gene_id (they were together in my data)
  #for example: ENSG00000234711|TUBB8P11 -> ENSG00000234711(in gene_id column) , TUBB8P11 (in gene_name column)
  df$gene_name <- str_sub(df$gene_names,17,-1)
  df$gene_id <- str_sub(df$gene_names, 1, 15)
  df$gene_names <- NULL
  #reordering columns (this should be checked!)
  combination = c(df_column_number+2, df_column_number+1, seq(df_column_number))
  df <- df[,combination]
  #setting gene_id to gene_names for the ones that are empty
  #for example some genes do not have gene names only gene id: ENSG00000288531, in this case this is set as gene_name
  #for loops over rows
  for(i in 1:nrow(df)) {      
    if (df[i ,2 ] ==""){
      df[i,2] <- df[i,1]
    }
  }
  return (df)
}
#making resSig **better**
resSig <-data_frame_me(resSig)

#getting the normalized&transformed matrix (same working scheme compared to resSig)
vsd_mat_df <- data_frame_me(vsd_mat)

#getting the normalized counts for ONLY significant genes (then we will create a heatmap)
resSig <- merge(resSig, vsd_mat_df, by =c("gene_id","gene_name"))

#exporting as xlsx
write_xlsx(resSig,"WT_vs_DGEZ_differential_expressed_genes.xlsx")

#we do not need the whole dataset, to use it in a graph, slicing is done
resSig_graph <- resSig[,c(2,8:length(colnames(resSig)))]
rownames(resSig_graph) <- resSig_graph$gene_name
resSig_graph$gene_name <- NULL

#setting the column names of this new dataframe
#colnames(resSig_graph) <- c("PATZ1_G1-2","PATZ1_G1-3","PATZ1_WT1", "PATZ1_WT2", "PATZ1_WT3")

pheatmap(resSig_graph, cluster_rows=TRUE, show_rownames=FALSE,scale = "row",
         cluster_cols=TRUE, annotation_col=select(metadata, condition),file="WT_vs_DGE_heatmap.tiff",width=10,height=10,
         main = "Wild Type vs DGEZ Differentially Expressed Genes Heatmap")

#getting the lab argument for volcano plot (it requires ALL gene names, same structure as previous two
result_shr_df <- data_frame_me(result_shr)
#getting all gene counts as well (just for export)
just_for_export <- merge(result_shr_df, vsd_mat_df, by =c("gene_id","gene_name"))
#exporting as xlsx
write_xlsx(just_for_export,"WT_vs_DGEZ_all_genes.xlsx")


jpeg(file="WT_vs_DGEZ_volcano.tiff",width = 6000, height= 6000,res = 600)

#Doing volcano plot
EnhancedVolcano(result_shr_df,
                lab = result_shr_df$gene_name,
                xlim = c(-10,10),
                pCutoff = 5e-02,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'padj',
                title = "Wild Type vs DGEZ Differentially Expressed Genes")

dev.off()
#comparing the results to patz2 (from another set of experiments) 
#VENN diagram will be used
library("readxl")
library(VennDiagram)

WT_vs_DGEZ <- read_excel("WT_vs_DGEZ_differential_expressed_genes.xlsx")
WT_vs_DGE<- read_excel("WT_vs_DGE_differential_expressed_genes.xlsx")
WT_vs_iT2puro <- read_excel("WT_vs_iT2puro_differential_expressed_genes.xlsx")


#getting only the gene_ids from both parts and setting to set1 and set2 variables

set1 <- WT_vs_DGEZ$gene_id
set2 <- WT_vs_DGE$gene_id 
set3 <- WT_vs_iT2puro$gene_id


#we will see the COMMON genes beforehand
library(dplyr)
DGEZ_unique_genes <- anti_join(WT_vs_DGEZ, WT_vs_iT2puro, by=c("gene_id","gene_name"))
DGE_unique_genes <- anti_join(WT_vs_DGE, WT_vs_iT2puro, by=c("gene_id","gene_name"))
DGE_u_DGEZ_unique_genes <- merge(DGEZ_unique_genes, DGE_unique_genes, by=c("gene_id","gene_name"))

write_xlsx(DGEZ_unique_genes,"DGEZ_unique_genes.xlsx")
write_xlsx(DGE_unique_genes,"DGE_unique_genes.xlsx")
write_xlsx(DGE_u_DGEZ_unique_genes,"DGE_u_DGEZ_unique_genes.xlsx")


#Now to the venn diagram (look further for venn.diagram manual)
#set1 and set2 are used, category.names can be changed and this automatically saves a file (filename argument)
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Wild Type vs DGEZ" , "Wild Type vs DGE", "Wild Type vs iT2puro"),
  filename="WT_vs_DGEZ_vs_DGE_vs_i2Tpuro_DEGS.png",
  output=TRUE,
  
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = c("darkseagreen","coral2","lavender"),
  
  # Numbers
  cex = 1.2,
  fontface = "bold",
  fontfamily = "sans",
  
  #set names
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 180),
  cat.fontfamily = "sans",
)

#getting the toP10 significant genes (order with padj values)
top20 <- resSig[order(resSig$padj), ]
top20 <- top20[1:20,c(2,8:12)]

#now to convert this dataframe into a key-value pairs with gather function
top20 <- gather(top20, key = "Sample Name", value= "nrm_counts", 2:6)

#inner joining metadata with the toP10 variable (inserting metadata information to toP10 df)
top20 <- inner_join(top20, rownames_to_column(metadata, var = "Sample Name"), by="Sample Name")

jpeg(file="WT_vs_DGEZ.tiff",width = 6000, height= 4000,res = 600)

#creating the expression plot
ggplot (top20) +
  geom_point(aes(x=gene_name, y= nrm_counts,color =condition)) +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Most Significant Differentially Expressed Genes for Wild Type vs DGEZ")+
  theme(axis.text.x = element_text(angle = 45, hjust =1))+
  theme(plot.title = element_text(hjust=0.5))

dev.off()
