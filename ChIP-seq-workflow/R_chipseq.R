if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("ChIPQC")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("ChIPseeker")
BiocManager::install("EnsDb.Hsapiens.v86")

install.packages("rlang")
install.packages("UpSetR")
install.packages("ggupset")
install.packages("ggimage")
install.packages( "Nozzle.R1", type="source" )

library(ChIPQC)
library(BiocParallel)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(UpSetR)
library(ggupset)
library(ggimage)
library(Nozzle.R1)

setwd("C:/Users/sofiap/R/chip_seq")

#creating the metadata (samples)
#the names of the variables MUST be like those:
samples <- read.csv('meta/metadata.txt')

#if you have only one sample use this instead
#samples <- ChIPQCsample("data/bams/kaiso.bam")

#creating chipobj, the indexes of the bam files referred in the bamReads should be given
register(SerialParam()) #parallel isn't working, so serial is done

#consensus and bCount,summits are from DiffBind
chipObj <- ChIPQC(samples,annotation = "hg38",consensus=TRUE, bCount=TRUE, summits=250)

QCmetrics(chipObj)

#creating a report

ChIPQCreport(chipObj, reportName="ChIP QC report: ZBTB24 Chip", reportFolder="results/zbtb24_chip")

#!!!!! ChIPseeker PART !!!!!

#samplefiles list is constructed, note the order of names

samplefiles <- list.files("data/peakcalls/idr", pattern = ".bed", full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("zbtb24_all","zbtb24_idr")

#you can extract the peaks to separate variables
peak_zbtb24 <- readPeakFile(samplefiles[[2]])
#peak_kaiso <- readPeakFile(samplefiles[[2]])

#coverage plotting
covplot(peak_zbtb24,weightCol = "V5",chrs=c("chr17", "chr18"))


#assign annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #if you used uscs as a reference use this
edb <- EnsDb.Hsapiens.v86 #if you used ensembl as a reference genome use this 

#annotation of peak
peakAnno <- annotatePeak(peak_zbtb24, tssRegion=c(-1000, 1000), TxDb=edb, annoDb="org.Hs.eg.db")


#annotation of peaks (if you want to do more than one at the same time)
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=edb, annoDb = "org.Hs.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

#plotting the annotation distribution
plotAnnoBar(peakAnno) #barplot

#plotAnnoPie(peakAnnoList[[1]]) #pieplot
plotAnnoPie(peakAnno)

#cool plot
upsetplot(peakAnno, vennpie=TRUE) #upsetplot

#plotting to TF binding loci relative to the transcription starting site
plotDistToTSS(peakAnno, title="Distribution of Transcription factor-binding loci relative to TSS")

#plotting average profiles
promoter <- getPromoters(TxDb=edb, upstream=1000, downstream=1000)
tagMatrixList <- lapply(samplefiles, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList[[2]], xlim=c(-1000, 1000),conf=0.95,resample=500, facet="row")

#peak profiles
#plotPeakProf2(tagMatrixList, upstream = rel(0.2), downstream = rel(0.2),
#              conf = 0.95, by = "gene", type = "body",
#              TxDb = edb, facet = "row", nbin = 800)


#you can extract the annotation dataframes one by one
zbtb24_anno <- data.frame(peakAnno@anno)
#kaiso_annot <- data.frame(peakAnnoList[["kaiso"]]@anno)


#lets get the common genes! (only from promoter sites)
zbtb24_anno <- subset(zbtb24_anno, annotation == "Promoter")
zbtb24_anno <- subset(zbtb24_anno, transcriptBiotype == "protein_coding")

#save your dataframe if wanted
write.csv(zbtb24_anno, "zbtb24_idr_anno_promoter_protein_coding.csv")


#kaiso_annot <- subset(kaiso_annot, annotation == "Promoter")

zbtb24_genes <- unique(zbtb24_anno$geneId)
bcl6_genes <- bcl6_genes[!is.na(bcl6_genes)]

kaiso_genes <- unique(kaiso_annot$ENTREZID)
kaiso_genes <- kaiso_genes[!is.na(kaiso_genes)]

#enrichment analysis KEGG

genes = lapply(peakAnno, function(i) as.data.frame(i)$ENTREZID)
names(genes) = sub("_", "\n", names(genes))

#or do this (they should be ENTREZID)
zbtb24_genes <- (unique(zbtb24_chip[["ENTREZID"]]))
patz1_genes <- na.omit((as.integer(patz1_deg[["entrez"]])))

gene_list <- list(patz1_genes,zbtb24_genes)
names(gene_list) <- c("PATZ1","ZBTB24")
compKEGG <- compareCluster(geneCluster   = gene_list,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.1,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
