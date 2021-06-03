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
wd <- "/Users/haohu/Research/AZ-RNASeq/ME34"
wd <- getwd()

list.files(wd)
#prepare sample table for DEseq analysis
ME34Files <- grep('HTseq-count_hisat2_ME34.*.out$', list.files(wd), value=TRUE)
ME34Files
ME34Names <- gsub("HTseq-count_hisat2_(.*)_[A|G|C|T]{6}.out", "\\1", ME34Files)
ME34Names
ME34Conditions <- gsub("HTseq-count_hisat2_(.*)_([A|L|U])([1-3])_[A|G|C|T]{6}.out", "\\1-\\2-\\3", ME34Files)
ME34Conditions
ME34Table <- data.frame(sampleName=ME34Names, 
                          fileName=ME34Files, 
                          condition=ME34Conditions)
ME34Table <- separate(ME34Table, condition, into = c("time","tissue", "batch"), sep = "-")
ME34Table$time <- factor(ME34Table$time, levels = c("ME34_young", "ME34_old"))
ME34Table$tissue <- factor(ME34Table$tissue, levels = c("L", "A", "U"))
ME34Table$batch <- factor(ME34Table$batch, levels = c(1, 2, 3))
#group the time and tissue variable to reduce factors
ME34Table$group <- factor(paste0(ME34Table$time, ME34Table$tissue))
ME34Table

####### input table into DESeq2 and DESeq2 analysis ######
ME34HTseq <- DESeqDataSetFromHTSeqCount(sampleTable = ME34Table, 
                                        directory = wd, 
                                        design = ~ batch + group)

####analysis with filtering sum of row counts>=10 and at least one count > 10

ME34HTseq1 <- ME34HTseq[rowSums(counts(ME34HTseq)) >=10 & (rowSums(counts(ME34HTseq) > 10) >= 1), ]
dds_ME34 <- DESeq(ME34HTseq1)
plotDispEsts(dds_ME34)

#get the normalized gene counts into a dataframe
counts_ME34 <- as.data.frame(counts(dds_ME34, normalized = TRUE))

#add a new column of gene names
counts_ME34$geneID <- rownames(counts_ME34)

#extracting transformed values
vst1 <- vst(dds_ME34, blind = FALSE)
rld <- rlog(dds_ME34, blind = FALSE)
rld
head(assay(rld))
meanSdPlot(assay(rld), ranks = FALSE)

ME34_rld <- as.data.frame(assay(rld))
ME34_rld$Sv_geneID <- rownames(ME34_rld)

write.table(ME34_rld, "ME34_rld.txt", sep = "\t", row.names = F)

#PCA graph
pcaData <- plotPCA(rld, intgroup = c("group", "batch"), returnData = TRUE)
pcaData
percentVar <- round(100*attr(pcaData, "percentVar"))
pdf("ME34_PCA.pdf", width = 5, height = 7)
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
pdf("ME34_sampledistance.pdf", width = 6, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
#pvclust
library(pvclust)
result <- pvclust(assay(rld))
pdf("ME34_pvclust.pdf", width = 8, height = 8)
plot(result)
dev.off()  

###### differential expression analysis ######
####run DE analysis
##ME34_young_A_vs_L
res_ME34_young_A_vs_L <- results(dds_ME34, 
                                 contrast = c("group", "ME34_youngA", "ME34_youngL")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_ME34_young_A_vs_L)
dim(res_ME34_young_A_vs_L)
mcols(res_ME34_young_A_vs_L)$description  #describe what each column does

#MA plot
pdf("ME34_young_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_ME34_young_A_vs_L, ylim = c(-8,8), main = 'ME34_young_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_ME34_young_A_vs_L <- identify(res_ME34_young_A_vs_L$baseMean, 
                                  res_ME34_young_A_vs_L$log2FoldChange)
rownames(res_ME34_young_A_vs_L)[idx_ME34_young_A_vs_L]

##ME34_old_A_vs_L
res_ME34_old_A_vs_L <- results(dds_ME34, 
                               contrast = c("group", "ME34_oldA", "ME34_oldL")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_ME34_old_A_vs_L)
dim(res_ME34_old_A_vs_L)


#MA plot
pdf("ME34_old_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_ME34_old_A_vs_L, ylim = c(-8,8), main = 'ME34_old_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_ME34_old_A_vs_L <- identify(res_ME34_old_A_vs_L$baseMean, 
                                res_ME34_old_A_vs_L$log2FoldChange)
rownames(res_ME34_old_A_vs_L)[idx_ME34_old_A_vs_L]

##ME34_young_A_vs_U
res_ME34_young_A_vs_U <- results(dds_ME34, 
                                 contrast = c("group", "ME34_youngA", "ME34_youngU")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_ME34_young_A_vs_U)
dim(res_ME34_young_A_vs_U)


#MA plot
pdf("ME34_young_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_ME34_young_A_vs_U, ylim = c(-8,8), main = 'ME34_young_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_ME34_young_A_vs_U <- identify(res_ME34_young_A_vs_U$baseMean, 
                                  res_ME34_young_A_vs_U$log2FoldChange)
rownames(res_ME34_young_A_vs_U)[idx_ME34_young_A_vs_U]

##ME34_old_A_vs_U
res_ME34_old_A_vs_U <- results(dds_ME34, 
                               contrast = c("group", "ME34_oldA", "ME34_oldU")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_ME34_old_A_vs_U)
dim(res_ME34_old_A_vs_U)


#MA plot
pdf("ME34_old_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_ME34_old_A_vs_U, ylim = c(-8,8), main = 'ME34_old_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_ME34_old_A_vs_U <- identify(res_ME34_old_A_vs_U$baseMean, 
                                res_ME34_old_A_vs_U$log2FoldChange)
rownames(res_ME34_old_A_vs_U)[idx_ME34_old_A_vs_U]

##ME34_old_vs_young_A
res_ME34_old_vs_young_A <- results(dds_ME34, 
                           contrast = c("group", "ME34_oldA", "ME34_youngA")) 
                           # lfcThreshold = 0.585,
                           # alpha = 0.05,
                           # altHypothesis = "greaterAbs")
summary(res_ME34_old_vs_young_A)
dim(res_ME34_old_vs_young_A)
# 
#MA plot
pdf("ME34_old_vs_young_A.pdf", width = 4, height = 4)
plotMA(res_ME34_old_vs_young_A, ylim = c(-12,12), main = 'ME34_old_vs_young_A')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()

####join the res_tables together
##join the DE table and gene count tables together
names(res_ME34_young_A_vs_L) <-paste0(names(res_ME34_young_A_vs_L),"_ME34youngAvsL")
names(res_ME34_young_A_vs_U) <-paste0(names(res_ME34_young_A_vs_U),"_ME34youngAvsU")
names(res_ME34_old_A_vs_L) <-paste0(names(res_ME34_old_A_vs_L),"_ME34oldAvsL")
names(res_ME34_old_A_vs_U) <-paste0(names(res_ME34_old_A_vs_U),"_ME34oldAvsU")
names(res_ME34_old_vs_young_A) <-paste0(names(res_ME34_old_vs_young_A),"_ME34youngvsoldA")

res_ME34 <- cbind(res_ME34_young_A_vs_L, res_ME34_young_A_vs_U,
                  res_ME34_old_A_vs_L, res_ME34_old_A_vs_U, 
                  res_ME34_old_vs_young_A, counts_ME34)
dim(res_ME34)
names(res_ME34)
names(res_ME34)[49] <- 'Sv_geneID'

##input S.viridis annotation
Sv_anno <- fread("/Users/haohu/Research/Genomes/Sviridis_311_v1.1.annotation_info.txt", 
                 header = T, sep = "\t", select = c(1:16))%>%tbl_df()
dim(Sv_anno)
names(Sv_anno)
Sv_anno <-Sv_anno %>%
  dplyr::select(geneID = locusName, Best_hit_arabi_name:rice_defline, Pfam:GO)
Sv_anno <-Sv_anno[!duplicated(Sv_anno$geneID),]
dim(Sv_anno)
names(Sv_anno) <- paste0("Sv_", names(Sv_anno))
head(Sv_anno)

##join the result table and annotation together
res_ME34_anno <- left_join(as.data.frame(res_ME34), Sv_anno, by = "Sv_geneID")
dim(res_ME34_anno)
write.table(res_ME34_anno, "res_ME34_anno_DEseq2_filter10.txt", row.names = F, sep = "\t")
