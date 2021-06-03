####### load the libraries ######
library(DESeq2)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(writexl)

###### import and format table for DESeq2 ######
wd <- "/Users/haohu/Research/AZ-RNASeq/Brachy"
wd <- getwd()

list.files(wd)
#prepare sample table for DEseq analysis
BrachyFiles <- grep('HTseq-count_hisat2_Brachy.*.out$', list.files(wd), value=TRUE)
BrachyFiles
BrachyNames <- gsub("HTseq-count_hisat2_(.*)_[A|G|C|T]{6}.out", "\\1", BrachyFiles)
BrachyNames
BrachyConditions <- gsub("HTseq-count_hisat2_(.*)_([A|L|U])([1-3])_[A|G|C|T]{6}.out", "\\1-\\2-\\3", BrachyFiles)
BrachyConditions
BrachyTable <- data.frame(sampleName=BrachyNames, 
                          fileName=BrachyFiles, 
                          condition=BrachyConditions)
BrachyTable <- separate(BrachyTable, condition, into = c("time","tissue", "batch"), sep = "-")
BrachyTable$time <- factor(BrachyTable$time, levels = c("Brachy_young", "Brachy_old"))
BrachyTable$tissue <- factor(BrachyTable$tissue, levels = c("L", "A", "U"))
BrachyTable$batch <- factor(BrachyTable$batch, levels = c(1, 2, 3))
#group the time and tissue variable to reduce factors
BrachyTable$group <- factor(paste0(BrachyTable$time, BrachyTable$tissue))
BrachyTable

####### input table into DESeq2 and DESeq2 analysis ######
BrachyHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = BrachyTable, 
                                        directory = wd, 
                                        design = ~ batch + group)

####analysis with filtering sum of row counts>=10 and at least one count > 10

BrachyHTseq1 <- BrachyHTseq[rowSums(counts(BrachyHTseq)) >=10 & (rowSums(counts(BrachyHTseq) > 10) >= 1), ]
dds_Brachy <- DESeq(BrachyHTseq1)
plotDispEsts(dds_Brachy)

#get the normalized gene counts into a dataframe
counts_Brachy <- as.data.frame(counts(dds_Brachy, normalized = TRUE))

#add a new column of gene names
counts_Brachy$geneID <- rownames(counts_Brachy)

#extracting transformed values
vst1 <- vst(dds_Brachy, blind = FALSE)
rld <- rlog(dds_Brachy, blind = FALSE)
rld
head(assay(rld))
meanSdPlot(assay(rld), ranks = FALSE)

Brachy_rld <- as.data.frame(assay(rld))
Brachy_rld$Bd_geneID <- rownames(Brachy_rld)

write.table(Brachy_rld, "Brachy_rld.txt", sep = "\t", row.names = F)

#PCA graph
pcaData <- plotPCA(rld, intgroup = c("group", "batch"), returnData = TRUE)
pcaData
percentVar <- round(100*attr(pcaData, "percentVar"))
pdf("Brachy_PCA.pdf", width = 5, height = 7)
ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = batch)) + geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()

#sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf("Brachy_sampledistance.pdf", width = 6, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
#pvclust
library(pvclust)
result <- pvclust(assay(rld))
pdf("Brachy_pvclust.pdf", width = 8, height = 8)
plot(result)
dev.off()  

###### differential expression analysis ######
####run DE analysis
##Brachy_young_A_vs_L
res_Brachy_young_A_vs_L <- results(dds_Brachy, 
                                 contrast = c("group", "Brachy_youngA", "Brachy_youngL")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_Brachy_young_A_vs_L)
dim(res_Brachy_young_A_vs_L)
mcols(res_Brachy_young_A_vs_L)$description  #describe what each column does

#MA plot
pdf("Brachy_young_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_Brachy_young_A_vs_L, ylim = c(-8,8), main = 'Brachy_young_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Brachy_young_A_vs_L <- identify(res_Brachy_young_A_vs_L$baseMean, 
                                  res_Brachy_young_A_vs_L$log2FoldChange)
rownames(res_Brachy_young_A_vs_L)[idx_Brachy_young_A_vs_L]

##Brachy_old_A_vs_L
res_Brachy_old_A_vs_L <- results(dds_Brachy, 
                               contrast = c("group", "Brachy_oldA", "Brachy_oldL")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_Brachy_old_A_vs_L)
dim(res_Brachy_old_A_vs_L)


#MA plot
pdf("Brachy_old_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_Brachy_old_A_vs_L, ylim = c(-8,8), main = 'Brachy_old_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Brachy_old_A_vs_L <- identify(res_Brachy_old_A_vs_L$baseMean, 
                                res_Brachy_old_A_vs_L$log2FoldChange)
rownames(res_Brachy_old_A_vs_L)[idx_Brachy_old_A_vs_L]

##Brachy_young_A_vs_U
res_Brachy_young_A_vs_U <- results(dds_Brachy, 
                                 contrast = c("group", "Brachy_youngA", "Brachy_youngU")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_Brachy_young_A_vs_U)
dim(res_Brachy_young_A_vs_U)


#MA plot
pdf("Brachy_young_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_Brachy_young_A_vs_U, ylim = c(-8,8), main = 'Brachy_young_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Brachy_young_A_vs_U <- identify(res_Brachy_young_A_vs_U$baseMean, 
                                  res_Brachy_young_A_vs_U$log2FoldChange)
rownames(res_Brachy_young_A_vs_U)[idx_Brachy_young_A_vs_U]

##Brachy_old_A_vs_U
res_Brachy_old_A_vs_U <- results(dds_Brachy, 
                               contrast = c("group", "Brachy_oldA", "Brachy_oldU")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_Brachy_old_A_vs_U)
dim(res_Brachy_old_A_vs_U)


#MA plot
pdf("Brachy_old_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_Brachy_old_A_vs_U, ylim = c(-8,8), main = 'Brachy_old_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Brachy_old_A_vs_U <- identify(res_Brachy_old_A_vs_U$baseMean, 
                                res_Brachy_old_A_vs_U$log2FoldChange)
rownames(res_Brachy_old_A_vs_U)[idx_Brachy_old_A_vs_U]

##Brachy_old_vs_young_A
res_Brachy_old_vs_young_A <- results(dds_Brachy, 
                           contrast = c("group", "Brachy_oldA", "Brachy_youngA")) 
                           # lfcThreshold = 0.585,
                           # alpha = 0.05,
                           # altHypothesis = "greaterAbs")
summary(res_Brachy_old_vs_young_A)
dim(res_Brachy_old_vs_young_A)
# 
#MA plot
pdf("Brachy_old_vs_young_A.pdf", width = 4, height = 4)
plotMA(res_Brachy_old_vs_young_A, ylim = c(-12,12), main = 'Brachy_old_vs_young_A')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()

####join the res_tables together
##join the DE table and gene count tables together
names(res_Brachy_young_A_vs_L) <-paste0(names(res_Brachy_young_A_vs_L),"_BrachyyoungAvsL")
names(res_Brachy_young_A_vs_U) <-paste0(names(res_Brachy_young_A_vs_U),"_BrachyyoungAvsU")
names(res_Brachy_old_A_vs_L) <-paste0(names(res_Brachy_old_A_vs_L),"_BrachyoldAvsL")
names(res_Brachy_old_A_vs_U) <-paste0(names(res_Brachy_old_A_vs_U),"_BrachyoldAvsU")
names(res_Brachy_old_vs_young_A) <-paste0(names(res_Brachy_old_vs_young_A),"_BrachyyoungvsoldA")

res_Brachy <- cbind(res_Brachy_young_A_vs_L, res_Brachy_young_A_vs_U,
                  res_Brachy_old_A_vs_L, res_Brachy_old_A_vs_U, 
                  res_Brachy_old_vs_young_A, counts_Brachy)
dim(res_Brachy)
names(res_Brachy)
names(res_Brachy)[49] <- 'Bd_geneID'

##input S.viridis annotation
Bd_anno <- fread("/Users/haohu/Research/Genomes/Bdistachyon_314_v3.1.annotation_info.txt", 
                 header = T, sep = "\t", select = c(1:16))%>%tbl_df()
dim(Bd_anno)
names(Bd_anno)
Bd_anno <-Bd_anno %>%
  dplyr::select(geneID = locusName, Best_hit_arabi_name:rice_defline, Pfam:GO)
Bd_anno <-Bd_anno[!duplicated(Bd_anno$geneID),]
dim(Bd_anno)
names(Bd_anno) <- paste0("Bd_", names(Bd_anno))
head(Bd_anno)

##join the result table and annotation together
res_Brachy_anno <- left_join(as.data.frame(res_Brachy), Bd_anno, by = "Bd_geneID")
dim(res_Brachy_anno)
write.table(res_Brachy_anno, "res_Brachy_anno_DEseq2_filter10.txt", row.names = F, sep = "\t")
