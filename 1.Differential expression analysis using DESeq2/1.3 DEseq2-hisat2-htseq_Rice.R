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
wd <- "/Users/haohu/Research/AZ-RNASeq/Rice"
wd <- getwd()

list.files(wd)
#prepare sample table for DEseq analysis
RiceFiles <- grep('HTseq-count_hisat2_Rice.*.out$', list.files(wd), value=TRUE)
RiceFiles
RiceNames <- gsub("HTseq-count_hisat2_(.*)_[A|G|C|T]{6}.out", "\\1", RiceFiles)
RiceNames
RiceConditions <- gsub("HTseq-count_hisat2_(.*)_([A|L|U])([1-3])_[A|G|C|T]{6}.out", "\\1-\\2-\\3", RiceFiles)
RiceConditions
RiceTable <- data.frame(sampleName=RiceNames, 
                          fileName=RiceFiles, 
                          condition=RiceConditions)
RiceTable <- separate(RiceTable, condition, into = c("time","tissue", "batch"), sep = "-")
RiceTable$time <- factor(RiceTable$time, levels = c("Rice_young", "Rice_old"))
RiceTable$tissue <- factor(RiceTable$tissue, levels = c("L", "A", "U"))
RiceTable$batch <- factor(RiceTable$batch, levels = c(1, 2, 3))
#group the time and tissue variable to reduce factors
RiceTable$group <- factor(paste0(RiceTable$time, RiceTable$tissue))
RiceTable

####### input table into DESeq2 and DESeq2 analysis ######
RiceHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = RiceTable, 
                                        directory = wd, 
                                        design = ~ batch + group)

####analysis with filtering sum of row counts>=10 and at least one count > 10

RiceHTseq1 <- RiceHTseq[rowSums(counts(RiceHTseq)) >=10 & (rowSums(counts(RiceHTseq) > 10) >= 1), ]
dds_Rice <- DESeq(RiceHTseq1)
plotDispEsts(dds_Rice)

#get the normalized gene counts into a dataframe
counts_Rice <- as.data.frame(counts(dds_Rice, normalized = TRUE))

#add a new column of gene names
counts_Rice$geneID <- rownames(counts_Rice)

#extracting transformed values
vst1 <- vst(dds_Rice, blind = FALSE)
rld <- rlog(dds_Rice, blind = FALSE)
rld
head(assay(rld))
meanSdPlot(assay(rld), ranks = FALSE)

Rice_rld <- as.data.frame(assay(rld))
Rice_rld$Os_geneID <- rownames(Rice_rld)

write.table(Rice_rld, "Rice_rld.txt", sep = "\t", row.names = F)

#PCA graph
pcaData <- plotPCA(rld, intgroup = c("group", "batch"), returnData = TRUE)
pcaData
percentVar <- round(100*attr(pcaData, "percentVar"))
pdf("Rice_PCA.pdf", width = 5, height = 7)
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
pdf("Rice_sampledistance.pdf", width = 6, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
#pvclust
library(pvclust)
result <- pvclust(assay(rld))
pdf("Rice_pvclust.pdf", width = 8, height = 8)
plot(result)
dev.off()  

###### differential expression analysis ######
####run DE analysis
##Rice_young_A_vs_L
res_Rice_young_A_vs_L <- results(dds_Rice, 
                                 contrast = c("group", "Rice_youngA", "Rice_youngL")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_Rice_young_A_vs_L)
dim(res_Rice_young_A_vs_L)
mcols(res_Rice_young_A_vs_L)$description  #describe what each column does

#MA plot
pdf("Rice_young_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_Rice_young_A_vs_L, ylim = c(-8,8), main = 'Rice_young_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Rice_young_A_vs_L <- identify(res_Rice_young_A_vs_L$baseMean, 
                                  res_Rice_young_A_vs_L$log2FoldChange)
rownames(res_Rice_young_A_vs_L)[idx_Rice_young_A_vs_L]

##Rice_old_A_vs_L
res_Rice_old_A_vs_L <- results(dds_Rice, 
                               contrast = c("group", "Rice_oldA", "Rice_oldL")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_Rice_old_A_vs_L)
dim(res_Rice_old_A_vs_L)


#MA plot
pdf("Rice_old_A_vs_L.pdf", width = 4, height = 4)
plotMA(res_Rice_old_A_vs_L, ylim = c(-8,8), main = 'Rice_old_A_vs_L')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Rice_old_A_vs_L <- identify(res_Rice_old_A_vs_L$baseMean, 
                                res_Rice_old_A_vs_L$log2FoldChange)
rownames(res_Rice_old_A_vs_L)[idx_Rice_old_A_vs_L]

##Rice_young_A_vs_U
res_Rice_young_A_vs_U <- results(dds_Rice, 
                                 contrast = c("group", "Rice_youngA", "Rice_youngU")) 
                                 # lfcThreshold = 0.585,
                                 # alpha = 0.05,
                                 # altHypothesis = "greaterAbs")
summary(res_Rice_young_A_vs_U)
dim(res_Rice_young_A_vs_U)


#MA plot
pdf("Rice_young_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_Rice_young_A_vs_U, ylim = c(-8,8), main = 'Rice_young_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Rice_young_A_vs_U <- identify(res_Rice_young_A_vs_U$baseMean, 
                                  res_Rice_young_A_vs_U$log2FoldChange)
rownames(res_Rice_young_A_vs_U)[idx_Rice_young_A_vs_U]

##Rice_old_A_vs_U
res_Rice_old_A_vs_U <- results(dds_Rice, 
                               contrast = c("group", "Rice_oldA", "Rice_oldU")) 
                               # lfcThreshold = 0.585,
                               # alpha = 0.05,
                               # altHypothesis = "greaterAbs")
summary(res_Rice_old_A_vs_U)
dim(res_Rice_old_A_vs_U)


#MA plot
pdf("Rice_old_A_vs_U.pdf", width = 4, height = 4)
plotMA(res_Rice_old_A_vs_U, ylim = c(-8,8), main = 'Rice_old_A_vs_U')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()
#identify the gene of interest on MA plot by clicking the dots on the graph
idx_Rice_old_A_vs_U <- identify(res_Rice_old_A_vs_U$baseMean, 
                                res_Rice_old_A_vs_U$log2FoldChange)
rownames(res_Rice_old_A_vs_U)[idx_Rice_old_A_vs_U]

##Rice_old_vs_young_A
res_Rice_old_vs_young_A <- results(dds_Rice, 
                           contrast = c("group", "Rice_oldA", "Rice_youngA")) 
                           # lfcThreshold = 0.585,
                           # alpha = 0.05,
                           # altHypothesis = "greaterAbs")
summary(res_Rice_old_vs_young_A)
dim(res_Rice_old_vs_young_A)
# 
#MA plot
pdf("Rice_old_vs_young_A.pdf", width = 4, height = 4)
plotMA(res_Rice_old_vs_young_A, ylim = c(-12,12), main = 'Rice_old_vs_young_A')
abline(h=c(-0.585,0.585), col="dodgerblue", lwd=2)
dev.off()

####join the res_tables together
##join the DE table and gene count tables together
names(res_Rice_young_A_vs_L) <-paste0(names(res_Rice_young_A_vs_L),"_RiceyoungAvsL")
names(res_Rice_young_A_vs_U) <-paste0(names(res_Rice_young_A_vs_U),"_RiceyoungAvsU")
names(res_Rice_old_A_vs_L) <-paste0(names(res_Rice_old_A_vs_L),"_RiceoldAvsL")
names(res_Rice_old_A_vs_U) <-paste0(names(res_Rice_old_A_vs_U),"_RiceoldAvsU")
names(res_Rice_old_vs_young_A) <-paste0(names(res_Rice_old_vs_young_A),"_RiceyoungvsoldA")

res_Rice <- cbind(res_Rice_young_A_vs_L, res_Rice_young_A_vs_U,
                  res_Rice_old_A_vs_L, res_Rice_old_A_vs_U, 
                  res_Rice_old_vs_young_A, counts_Rice)
dim(res_Rice)
names(res_Rice)
names(res_Rice)[49] <- 'Os_geneID'

##input S.viridis annotation
Os_anno <- fread("/Users/haohu/Research/Genomes/anno_Rice_R498.txt", 
                 header = T, sep = "\t", select = c(1:16))%>%tbl_df()
dim(Os_anno)
names(Os_anno)
Os_anno <-Os_anno %>%
  dplyr::select(geneID = locusName, Best_hit_arabi_name:rice_defline, Pfam:GO)
Os_anno <-Os_anno[!duplicated(Os_anno$geneID),]
dim(Os_anno)
names(Os_anno) <- paste0("Os_", names(Os_anno))
head(Os_anno)

##join the result table and annotation together
res_Rice_anno <- left_join(as.data.frame(res_Rice), Os_anno, by = "Os_geneID")
dim(res_Rice_anno)
write.table(res_Rice_anno, "res_Rice_anno_DEseq2_filter10.txt", row.names = F, sep = "\t")
