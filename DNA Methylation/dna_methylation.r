version
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("RColorBrewer", "pheatmap","tidyverse","minfi","minfiData","remotes","shinyMethyl","shinyMethylData","limma"))
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
BiocManager::install("devtools","liftOver")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
remotes::install_github("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
BiocManager::install("liftOver")
BiocManager::install("DMRcate")
BiocManager::install("missMethyl")
install.packages("gt")

#loading those libraries
library(pheatmap)
library(gt)
library(tidyverse)
library(minfi)
library(limma)
library(missMethyl)
library(minfiData)
library(DMRcate)
library(shinyMethyl)
library(shinyMethylData)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(dplyr)
library(readr)
library(ggrepel)

library(devtools)
install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

#Setting the working directory
setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/EPICv2_samples")

#load example data
baseDir <- system.file("extdata", package = "minfiData")

#2 folders for each slide and a csv file samplesheet
list.files(baseDir)

#the folder contains green and red idat files, each slide here has 3 arrays so 6 files
list.files(file.path(baseDir, "5723646052"))

#you can read the directory that returns a dataframe
patients <- read.metharray.sheet(baseDir)

#now the all info is in the RG set
RGset <- read.metharray.exp(targets = patients)
#just to look at it as dataframe (extracts phenotype data)
pd <- pData(RGset)
pd[,1:4]

#you can also just look at the sample sheet, skip=7 skips first 7 rows
#sentrix_id is slide number and sentrix_position is the array number
patients2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), 
                     stringsAsFactors = FALSE, skip=7)
#but there is no basename column, that is why we can manually add it
patients2$Basename <- file.path(baseDir, patients2$Sentrix_ID, 
                               paste0(patients2$Sentrix_ID, 
                                      patients2$Sentrix_Position))
#you can check the annotations used 
annotation(RGset)
                     
#quality control!
qcReport(RGset, sampNames= pd$Sample_Name,
         sampGroups = pd$Sample_Group,pdf="QC_report_sample.pdf")

#you can also use shinyMethly package for further QC
#first create shinymethlyset from your RGset
ShinyMethylSet <- shinySummarize(RGset)

#then you can start the application
runShinyMethyl(ShinyMethylSet)

#additionally MDS plot can be drawn (uses Euclidian distances)
names <- pd$Sample_Name
groups <- pd$Sample_Group
mdsPlot(RGset, sampNames=names, sampGroups=groups,pch=1)

#density plot can be drawn independently
densityPlot(RGset, sampGroups=groups)

#density bean plot can be drawn independently
par(mar=c(5,6,4,2)) #for the alignment of the next figure
densityBeanPlot(RGset, sampNames=names, sampGroups=groups)

#now we can do filtering to get rid of non-specifics
detP <- detectionP(RGset) #detects the p values 
dim(detP)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
#filtering out the ones with p-value smaller than 0.01     
absent<- apply(detP, 1, function(x) sum(x>0.01)) #if none of the values are bigger than 0.01 it will return 0  
RGset_filtered <- RGset[absent==0,]
dim(RGset)
dim(RGset_filtered)

#we can d preprocessing here, after filtering!
ssNoob_filtered <- preprocessNoob(RGset_filtered,dyeCorr = TRUE,dyeMethod="single")

#getting methylated and unmethylated sites
methylated_sites <- getMeth(ssNoob_filtered)
unmethylated_sites <- getUnmeth(ssNoob_filtered)

#beta and M (logit of beta) values of the sites
beta_values <- getBeta(ssNoob_filtered)
m_values <- getM(ssNoob_filtered)

#we can plot betas by type(infiniium), only one column at a time can be plotted
plotBetasByType(ssNoob_filtered[,1], main = "R02C02")

#a quality control on methylation object
qc <- getQC(ssNoob_filtered)
ssNoob_filtered <- addQC(ssNoob_filtered, qc = qc)
plotQC(qc)

#if you want to look at specific genomic positions you can plot singleCpGs
cpgs <- c("cg00050873", "cg00212031", "cg26684946", "cg00128718") #examples
par(mfrow=c(2,2)) #to align the plot
par(mar=c(4,5,1,2)) #for the alignment of the next figure
plotCpg(ssNoob_filtered, cpg=cpgs, pheno=groups, type="categorical")

#quality control on methylation object, it does a series of controls and can remove outliers
#if you are going to do this, there is no need to do getQC, getSex or fixMethOutliers
qc_from_minfiQC <- minfiQC(ssNoob_filtered, fixOutliers = TRUE, verbose =TRUE)
plotQC(ssNoob_filtered)

#to find differentially methylated positions from methyl object (F-test for categorical)
dmp <- dmpFinder(ssNoob_filtered, pheno=groups, type="categorical")
sum(dmp$qval < 0.05, na.rm=TRUE) #how many sites are differentially methylated

#you can get the genomic annotation
annot <- as.data.frame(getAnnotation(ssNoob_filtered, dropNonMapping = TRUE))
#only get the annotations for the ones that are significant
significant_dmp <-as.data.frame(annot[rownames(annot) %in% rownames(dmp)[dmp$qval<0.05],])
#to look at a specific gene
CNIH3_gene <- significant_dmp[grep("CNIH3", significant_dmp$UCSC_RefGene_Name), ]#to plot the associated loci
cpgs <- rownames(CNIH3_gene)
par(mfrow=c(3,3)) #to align the plot
par(mar=c(2,2,2,2)) #for the alignment of the next figure
plotCpg(ssNoob_filtered, cpg=cpgs, pheno=groups, type="categorical")

#now we can map to genome to to obtain GenomicMethyl object
genomic_methly_set <- mapToGenome(ssNoob_filtered)

#we can obtain GenomicRatio Object from GenomicMetyl object (what argument for keeping beta and M values)
genomic_ratio_set <- ratioConvert(genomic_methly_set, what = "both")

##you can remove the the probes with SNPs 
genomic_ratio_set_filtered <- dropLociWithSnps(genomic_ratio_set)

#you can remove X and Y chromosomes from the data
#first get the annotation of the genome
#then pick the X and Y chromosomes
#keep <- !(featureNames(genomic_ratio_set_filtered) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
#table(keep)
#genomic_ratio_set_filtered <- genomic_ratio_set_filtered[keep,]

#you can calculate the beta values or the m values after this as well! 
#you can do the annotations here as well 
#annot <- getAnnotation(genomic_ratio_set,dropNonMapping = TRUE)

#------------------------------#

setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_merged")
baseDir <- "/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_merged"

list.files(baseDir)
patients<- read_csv("SampleSheet.csv", col_types=cols(.default = "c")) #it reads as a tibble, change the datastructure to df
patients$Basename <- file.path(baseDir, patients$slide, 
                               paste0(patients$slide, "_",
                                      patients$array))

RGset_450k <- read.metharray.exp(targets = patients)
pd <- pData(RGset)
names <- pd$Patient_Name
groups <- pd$Patient_Group
pd

patients$ID <- paste(patients$Final_Group,patients$Patient_Name,sep=".")
sampleNames(RGset_450k) <- patients$ID
RGset_450k

#if we want to merge
RGset <- combineArrays(RGset_EPICv1, RGset_450k, outType = "IlluminaHumanMethylation450k")
#we have to add annotations to the RGset object, both manifest(the chip information) and annotation(the genome version)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGset)["annotation"] = "20a1.hg38"

#For EPICv1 and hg38(this is not correct probably, maybe do a lift over from 19 to 28)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPIC"
#RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

#for 450k methylation this is done automatically?

annEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#the loci names are different than the epicv1 names but there is also a column in the annotation
#file that is called EPICv1_Loci, maybe that can be used to merge EPICv2&v1 samples?
manifestEPICv2 = getManifest(IlluminaHumanMethylationEPICv2manifest)

densityPlot(RGset, sampGroups=groups,pal = brewer.pal(8, "Dark2"),legend=FALSE)
legend("top", legend = levels(factor(groups)),
       text.col=brewer.pal(8,"Dark2"))
par(mar=c(3,3,3,3)) #for the alignment of the next figure

#I don't think this is working well ? The colors are always wrong
densityBeanPlot(RGset, sampGroups=patients$Patient_Group,sampNames=rev(names),pal = brewer.pal(3, "Dark2"))


#In this qc report, the bean plot is also wrong
qcReport(RGset, sampNames= patients$ID, sampGroups = patients$Patient_Group,
         pdf="QC_report_sample_EPICv1.pdf")


ShinyMethylSet <- shinySummarize(RGset)
runShinyMethyl(ShinyMethylSet)

detP <- detectionP(RGset)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?

#you can make a bar plot to look at the quality of the samples
pal <- brewer.pal(8,"Dark2")
pal <- brewer.pal(8,"Pastel2")
barplot(colMeans(detP), col=pal[factor(groups)], las=2,
        cex.names=1, cex.axis = 0.8,
        ylab="Mean detection p-values",names.arg=names,ylim=c(0,0.01)) #can add ylim here
abline(h=0.01,col="red") #this line is for 0.01 value
legend("top", legend=levels(factor(groups)), fill=pal,
       bty="n")
dev.off()

#getting rid of the low quality ones, if they exist in any of the samples
absent<- apply(detP, 1, function(x) sum(x>0.01))
RGset_filtered <- RGset[absent==0,]

ssNoob_filtered <- preprocessNoob(RGset_filtered,dyeCorr = TRUE,dyeMethod="single")

#after normalization we can look at the density plot
densityPlot(getBeta(ssNoob_filtered), sampGroups=groups,main="Density Plot of Beta Values", legend=FALSE)
legend("top", legend = levels(factor(groups)),
       text.col=brewer.pal(8,"Dark2"), bty="n",text.font=2)
dev.off()

qc <- getQC(ssNoob_filtered)
ssNoob_filtered <- addQC(ssNoob_filtered, qc = qc)
plotQC(qc)

methylated_sites <- getMeth(ssNoob_filtered)
unmethylated_sites <- getUnmeth(ssNoob_filtered)

beta_values <- getBeta(ssNoob_filtered)
m_values_new <- getM(ssNoob_filtered)

#maybe we can get back to here by saving the csv file.
write.table(m_values, file="m_values_full.csv", sep=",", row.names=TRUE)
m_values <- read.csv("m_values_combined.csv") #here use the read.csv function not read_csv..

#you can remove X and Y chromosomes from the data
#first get the annotation of the genome
#then pick the X and Y chromosomes
keep <- !(rownames(final_m_values) %in% ann450k_hg38$probeID[ann450k_hg38$CpG_chrm %in% c("chrX","chrY")])
final_m_values_filtered  <- final_m_values[keep,]

#examining the columns
head(m_values)
rownames(m_values)[grepl("cg25384689", rownames(m_values))]
as.data.frame(m_values)[rownames(as.data.frame(m_values)) == "cg25384689_TC22", ]

#lets look at PCA, plotMDS is from limma package, mdsPlot is from minfi, the labels is optional
plotMDS(m_values_filtered, top=10000, gene.selection="common", labels=patients$Patient_Name,cex= 0.8,
        col=pal[factor(patients$sex)])
legend("bottom", legend=levels(factor(patients$sex)), text.col=pal,
       bg="white", text.font=2,cex=0.8,bty="n")
dev.off()

#examining the higher dimensions of variance
par(mfrow=c(1,3))

plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(1,3))
legend("bottom",legend = levels(factor(patients$location)), text.col = pal,
       bg = "white", text.font = 2, bty = "n")
plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(2,3))
legend("bottom", legend=levels(factor(patients$location)), text.col=pal,
       bg="white", text.font=2,bty="n")
plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(3,4))
legend("bottomleft", legend=levels(factor(patients$location)), text.col=pal,
       bg="white", text.font=2,bty="n")

dev.off()

#hierarchical clustering heat map, in default uses euclidian distance
m_values_fixed <- m_values_filtered_wo_outlier
colnames(m_values_fixed) <- sub(".*\\.", "", colnames(m_values_fixed)) #getting rid of OA,RA, etc.
vsd_cor <- cor(m_values_fixed)
vsd_cor <- as.data.frame(vsd_cor)
patients_for_pheatmap <- as.data.frame(patients_wo_outlier)[,c(6,8),drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients_wo_outlier$Patient_Name
pheatmap(vsd_cor,annotation_col= patients_for_pheatmap,
         main = "Hierarchical Heatmap of Patient Samples", filename ="after_x_y_removal/pheatmap_corr_wo_outlier.jpeg",width= 10, height = 8)
graphics.off()


cpgs <- c("cg00050873") #we should add _BC21 etc for the EPIC v2... or else it doesn't work
par(mfrow=c(2,2)) #how many plots in one figure
par(mar=c(4,5,1,2)) #margins
plotCpg(ssNoob_filtered, cpg=cpgs, pheno=groups, type="categorical")

#to find differentially methylated positions from methyl object (F-test for categorical)
#we cant specify the specific comparisons, it compares all three of them, that is why limma is better
dmp <- dmpFinder(ssNoob_filtered, pheno=groups, type="categorical") 
sum(dmp$qval < 0.05, na.rm=TRUE) #how many sites are differentially methylated (0?)

# finding differentially methylated positions in limma
# this is the factor of interest
cellType <- factor(patients$Final_Group)
#this is the individual effect (we don't have biological replicates.. is this needed?)
individual <- factor(patients$Patient_Name)
#creating a design matrix
design <- model.matrix(~0+cellType, data=patients)
colnames(design) <- c(levels(cellType))
#rownames(design) <- names #this crashes the system??
fit <- lmFit(final_m_values_filtered, design)
#be careful about the x vs y, the controls are x most of the time
#but the logFC will be x to y difference not y to x, maybe we can do it the opposite way? 
contMatrix <- makeContrasts(VeRA-Healthy,RA-Healthy,OA-Healthy,RA-VeRA,RA-OA, levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#-----------------------------------------#

#then getting the annotations
#if you want to do for 450k and EPICv1 maybe we can just get the ones from EPICv2 annotation because it is hg38
cg22779972 <- as.data.frame(annEPICv2)
non_unique_names <- annEPICv2_df$Methyl450_Loci[duplicated(annEPICv2_df$Methyl450_Loci)]
annEPICv2_df_uni <- annEPICv2_df[!annEPICv2_df$Methyl450_Loci %in% non_unique_names, ]
rownames(annEPICv2_df_uni) <- annEPICv2_df_uni$Methyl450_Loci
#half of it is not even unique - there are duplicates and also some portion does not even have a Methyl450 counterpart

get_annots <- function(df,annot) {
  # Set the remaining unique values as row names
  df <- as.data.frame(merge(annot, df, by = "row.names",all.y=TRUE))
  df <- df[order(df$adj.P.Val), ]
  #writing the table as a csv file
  rownames(df) <- df$Row.names
  df <- df[,-1]
  return(df)
}

#getting the top genes (there is a problem with genelist, manifest file doesn't match one hundred percent)
DMPs_resolv_vs_RA <- topTable(fit2,  num=10, coef=1 ,sort.by="P")
DMPs_resolv_vs_RA <- get_annots(DMPs_resolv_vs_RA,annEPICv2_df_uni)
write.table(DMPs_resolv_vs_RA, file="DMPs_resolv_vs_RA.csv",sep=",",row.names = FALSE)

DMPs_RA_vs_estRA <- topTable(fit2, num=10,coef=6, sort.by="P")
DMPs_RA_vs_estRA <- get_annots(DMPs_RA_vs_estRA,annEPICv2_df_uni)

#---------------------------------------#

#or maybe just get the hg38(the lifted over version from hg18) for 450k :)
ann450k_hg38 <- read.delim("~/Desktop/PhD/2024 First Term/DNA Methylation/450k_samples/HM450.hg38.manifest.gencode.v36.tsv", header = TRUE)
row.names(ann450k_hg38) <- ann450k_hg38$probeID
annEPIC_hg38 <- read.delim("EPIC.hg38.manifest.gencode.v36.tsv", header = TRUE)
row.names(annEPIC_hg38) <- annEPIC_hg38$probeID

DMPs_Healthy_vs_VeRA <- topTable(fit2, num=Inf,coef=1, sort.by="P")
DMPs_Healthy_vs_VeRA <- get_annots(DMPs_Healthy_vs_VeRA,ann450k_hg38)
write.table(DMPs_Healthy_vs_VeRA, file="after_x_y_removal/DMPs_Healthy_vs_VeRA.csv",sep=",",row.names = FALSE)

DMPs_Healthy_vs_RA <- topTable(fit2, num=Inf,coef=2, sort.by="P")
DMPs_Healthy_vs_RA <- get_annots(DMPs_Healthy_vs_RA,ann450k_hg38)
write.table(DMPs_Healthy_vs_RA, file="after_x_y_removal/DMPs_Healthy_vs_RA.csv",sep=",",row.names = FALSE)

DMPs_Healthy_vs_OA <- topTable(fit2, num = Inf, coef=3, sort.by="P")
DMPs_Healthy_vs_OA <- get_annots(DMPs_Healthy_vs_OA,ann450k_hg38)
write.table(DMPs_Healthy_vs_OA, file="after_x_y_removal/DMPs_Healthy_vs_OA.csv",sep=",",row.names = FALSE)

DMPs_VeRA_vs_RA <- topTable(fit2, num = Inf, coef=4, sort.by="P")
DMPs_VeRA_vs_RA <- get_annots(DMPs_VeRA_vs_RA,ann450k_hg38)
write.table(DMPs_VeRA_vs_RA, file="after_x_y_removal/DMPs_VeRA_vs_RA.csv",sep=",",row.names = FALSE)

#---------------------------------------#

#look at the significant loci with graphs (we can find loci associated with specific genes)
par(mfrow=c(2,2))
sapply(DMPs_Healthy_vs_RA$probeID[1:4], function(cpg){
  plotCpg(m_values, cpg=cpg, pheno=patients$Final_Group, ylab = "Beta values")
})

dev.off()

genomic_methly_set <- mapToGenome(ssNoob_filtered)
genomic_ratio_set <- ratioConvert(genomic_methly_set, what = "both")

##you can remove the the probes with SNPs 
#genomic_ratio_set_filtered <- dropLociWithSnps(genomic_ratio_set)

# Create Venn diagram
set_1 <- rownames(DMPs_Healthy_vs_OA_sig)
set_2 <- rownames(DMPs_Healthy_vs_RA_sig)
venn.plot <- venn.diagram(
  x = list(set1 = set_1, set2 = set_2),
  category.names = c("Healthy vs OA", "Healthy vs RA"), filename=NULL,
    
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = c("#FFCBCB","#3C5B6F"),        cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",        cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
)

# Plot Venn diagram
grid.draw(venn.plot)


#heatmap for differentially methylated positions -> the common and significant ones
DMPs_Healthy_vs_RA_sig <- subset(DMPs_Healthy_vs_RA, DMPs_Healthy_vs_RA$adj.P.Val <= 0.05 & abs(DMPs_Healthy_vs_RA$logFC) >=1)
DMPs_Healthy_vs_OA_sig <- subset(DMPs_Healthy_vs_OA, DMPs_Healthy_vs_OA$adj.P.Val <= 0.05 & abs(DMPs_Healthy_vs_OA$logFC) >=1)
DMPs_Healthy_vs_VeRA_sig <- subset(DMPs_Healthy_vs_VeRA, DMPs_Healthy_vs_VeRA$adj.P.Val <= 0.05)
#only get the all three common loci
common_sig_loci <- Reduce(intersect, list(rownames(DMPs_Healthy_vs_RA_sig), rownames(DMPs_Healthy_vs_OA_sig)))
m_values_only_common_genes <- final_m_values_filtered[common_sig_loci,]

#after this be careful of the 43th and 44th position, did you get rid of the outlier or not?
#you can also look at the unique genes
unique_to_RA <- rownames(DMPs_Healthy_vs_RA_sig)[!(rownames(DMPs_Healthy_vs_RA_sig) %in% common_sig_loci)]

#or you can just look at the common genes
m_values_fixed <- m_values_only_common_genes #or unique genes
#m_values_fixed <- m_values_fixed[,c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42)] #order according to patient type
m_values_fixed <- m_values_fixed[,c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43)] #order according to joint location
colnames(m_values_fixed) <- sub(".*\\.", "", colnames(m_values_fixed)) #getting rid of OA,RA, etc.
patients_for_pheatmap <- as.data.frame(patients)[,c(8,6,4),drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients$Patient_Name
#patients_for_pheatmap <- patients_for_pheatmap[c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42),,drop=FALSE]
patients_for_pheatmap <- patients_for_pheatmap[c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43),,drop=FALSE]

pheatmap(m_values_fixed, cluster_rows=TRUE, show_rownames=FALSE,scale = "row",
         cluster_cols=TRUE, annotation_col=patients_for_pheatmap,file="after_x_y_removal/dmrs_heatmap_unique_to_RA_patient_location_rows_clustered.jpg",width=10,height=10,
         main = "Unique To RA Differentially Methylated Locations")

#only get the loci that are related to protein coding
common_sig_loci_protein <- Reduce(intersect, list(rownames(DMPs_Healthy_vs_OA[grepl("protein_coding",DMPs_Healthy_vs_OA$transcriptTypes),]), rownames(DMPs_Healthy_vs_RA[grepl("protein_coding",DMPs_Healthy_vs_RA$transcriptTypes),]), rownames(DMPs_Healthy_vs_VeRA[grepl("protein_coding",DMPs_Healthy_vs_VeRA$transcriptTypes),])))
m_values_only_common_genes_protein <- m_values_filtered[common_sig_loci_protein,]
m_values_fixed <- m_values_only_common_genes_protein
colnames(m_values_fixed) <- sub(".*\\.", "", colnames(m_values_fixed)) #getting rid of OA,RA, etc.
patients_for_pheatmap <- as.data.frame(patients)[,8,drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients$Patient_Name
pheatmap(m_values_fixed, cluster_rows=TRUE, show_rownames=FALSE,scale = "row",
         cluster_cols=TRUE, annotation_col=patients_for_pheatmap,file="after_x_y_removal/dmrs_only_proteins_heatmap.jpg",width=10,height=10,
         main = "Common Differentially Methylated Locations")

#volcano plot for comparison of healthy vs oa and healthy vs ra
combined_results <- merge(DMPs_Healthy_vs_OA_sig, DMPs_Healthy_vs_RA_sig, by = "probeID", suffixes = c("_OA", "_RA"), all = TRUE)
unique_to_OA <- combined_results[is.na(combined_results$logFC_RA), ]
unique_to_RA <- combined_results[is.na(combined_results$logFC_OA), ]
common <- combined_results[!is.na(combined_results$logFC_OA) & !is.na(combined_results$logFC_RA), ]

#create a coloumn describing uniqueness of significant loci
DMPs_Healthy_vs_RA$unique_col <- "N/A"
DMPs_Healthy_vs_RA$unique_col[DMPs_Healthy_vs_RA$probeID %in% unique_to_RA$probeID] <- "Unique to RA"
DMPs_Healthy_vs_RA$unique_col[DMPs_Healthy_vs_RA$probeID %in% common$probeID] <- "Common in RA & OA"

#create a column for gene labels
subset_df <- DMPs_Healthy_vs_RA[(DMPs_Healthy_vs_RA$unique_col == "Unique to RA" & grepl("protein_coding", DMPs_Healthy_vs_RA$transcriptTypes)),]
ordered_subset <- subset_df[order(subset_df$logFC), "probeID"]
ordered_subset_reverse <- subset_df[order(subset_df$logFC,decreasing=TRUE), "probeID"]
condition <- c(head(ordered_subset, 15),head(ordered_subset_reverse,15))
DMPs_Healthy_vs_RA$delabel <- ifelse(DMPs_Healthy_vs_RA$probeID %in% condition, DMPs_Healthy_vs_RA$genesUniq, NA)

# Plot the volcano plot
theme_set(theme_classic(base_size = 11) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            )
)
ggplot(data = DMPs_Healthy_vs_RA, aes(x = logFC, y = -log10(adj.P.Val), col= unique_col,label=delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size =1) +
  scale_color_manual(values = c("#E6A4B4", "gray", "#51829B"), 
                     labels = c("Common in RA & OA", "Not Significant", "Unique to RA")) +
  labs(color = 'Unique or Not') +
  geom_text_repel(max.overlaps = Inf) # To show all labels 
  


#Finding Differentially Methylated Regions!!DMRcate
#We will do lift-over from hg19 to hg38. I wont use lifted-over external hg38 because the CpGs are *possibly* not updated
library(liftOver)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#we extract the genome locations from the hg19
Human450klocs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations 
#Granges object is created so that we can do DMRcate 
Human450kRhg19 <- GRanges(paste(Human450klocs$chr, Human450klocs$pos, sep=":"))
names(Human450kRhg19) <- rownames(Human450klocs)

genome(Human450kRhg19) <- "hg19"
seqlevelsStyle(Human450kRhg19) <- "UCSC"

#getting the liftOver chain
ch <- import.chain("~/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/hg19ToHg38.over.chain")
Human450kRhg38 <- unlist(liftOver(Human450kRhg19, ch)) #this performs lifting over

#removing x and y from again because somehow some probes are retained
data.noXY <- rmSNPandCH(as.matrix(final_m_values_filtered), rmXY=TRUE)
#creating myAnnotation
myAnnotation <- cpg.annotate(object = data.noXY, datatype = "array", what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "VeRA - Healthy", arraytype = "450K",fdr=0.05)

#lifting over the myAnnotation object
myAnnotation.hg19 <- myAnnotation@ranges
retain <- names(myAnnotation.hg19) %in% names(Human450kRhg38)
myAnnotation.hg19 <- myAnnotation.hg19[retain]
myAnnotation.hg38 <- Human450kRhg38[names(myAnnotation.hg19)]
values(myAnnotation.hg38) <- values(myAnnotation.hg19)
#now we have our final annotation file!
myAnnotation <- new("CpGannotated", ranges=myAnnotation.hg38)

#DMRcate analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg38")
write.table(as.data.frame(results.ranges), file="after_x_y_removal/results.ranges_healthy_vs_VeRA.csv", sep=",", row.names=FALSE)
library(rtracklayer)
export(results.ranges,"results_ranges_healthy_vs_VeRA.bed",format="bed",index=FALSE)

# set up the grouping variables and colours
groups <- pal[1:length(unique(patients$Final_Group))]
names(groups) <- levels(factor(patients$Final_Group))
cols <- groups[as.character(factor(patients$Final_Group))]
# draw the plot for the top DMR
par(mfrow=c(1,1))
#works VERY slow
DMR.plot(ranges=results.ranges, dmr=1, CpGs=as.matrix(m_values_filtered), phen.col=cols, what = "M",
         arraytype = "450K",genome="hg38")

#b_values_are extracted as a bed file with 4 columns
#chr, chr start, chr end, value
#they are sorted, overlapping regions are averaged and returned
bed_file <- function(b_values_filtered,string){
df <- b_values_filtered[,grepl(string,colnames(b_values_filtered))]
df$avg <- rowMeans(df)*100
df <- merge(df, ann450k_hg38, by = "row.names", all.x = TRUE)
df <- df[,c("CpG_chrm","CpG_beg","CpG_end","avg")]
#b_values_healthy_vs_ra$name <- "."
#b_values_healthy_vs_ra$strand <- "."
#b_values_healthy_vs_ra <- b_values_healthy_vs_ra[, c(1,2,3,5,4,6)]
df <- na.omit(df)
df <- df[order(df$CpG_chrm, df$CpG_beg), ]
df <- aggregate(avg ~ CpG_chrm + CpG_beg + CpG_end, data = df, FUN = mean)
df <- df[order(df$CpG_chrm, df$CpG_beg), ] #should be sorted again after aggregation
return(df)
}
b_values_OA <- bed_file(b_values_filtered, "\\bOA\\.") #"\\bRA\\."
write.table(b_values_OA, file="after_x_y_removal/b_values_from_OA.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)

#GO Analysis
sigCpGs <- DMPs_Healthy_vs_RA$probeID[DMPs_Healthy_vs_RA$adj.P.Val<0.05]
all_loci <- DMPs_Healthy_vs_RA$probeID
#you can do this with a DMRs as well (results.ranges)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all_loci, plot.bias=TRUE)
top_pathways <- gst[order(gst$FDR), ]
top_pathways <- top_pathways[top_pathways$FDR < 0.05 & top_pathways$ONTOLOGY == "BP" ,]
top_pathways_second_col <- top_pathways[1:10,2]
top_pathways_rownames <- row.names(top_pathways)
new_df <- data.frame(rownames = top_pathways_rownames, first_column = top_pathways_second_col)
colnames(new_df) <- c("GO_Number","TERM")

new_df %>%gt() %>%
tab_header(title = "TOP 10 GO Biological Processes Pathways") %>%
tab_style(style = list(cell_fill(color = "#b2f7ef"),
cell_text(weight = "bold")),
locations = cells_body(columns = GO_Number))%>%
tab_style(style = list(cell_fill(color = "#ffefb5"),
cell_text(weight = "bold")), 
locations = cells_body(columns = TERM))
