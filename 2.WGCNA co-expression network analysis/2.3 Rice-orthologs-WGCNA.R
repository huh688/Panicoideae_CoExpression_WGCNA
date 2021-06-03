# Load the WGCNA package
library(WGCNA);

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()


###### Step 1 ###### Data input and cleaning
# Read in the data set
orthologs_expressed = read.csv("/Users/hao/Downloads/WGCNA/orthologs_expressed_counts.csv", header = TRUE);

# Take a quick look at what is in the data set:
dim(orthologs_expressed);
names(orthologs_expressed);

# Data exrtaction (normalized counts only) and transformation
datExpr = as.data.frame(t(orthologs_expressed[, c(36:38,33:35,39:41,27:29,24:26,30:32)]));
names(datExpr) = orthologs_expressed$Os_geneID;
dim(datExpr) 
rownames(datExpr)
colnames(datExpr)

#Check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

#Cluster the samples to identify any obvious outliers, and plot the  tree
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9) # Adjust dimensions if the window is too large or too small.

pdf(file = "Os_sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

###### Step 2 ###### Network construction and module detection
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
#pdf("soft threshold.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

# Check Scale-free topology
softPower = 16
k <- softConnectivity(datExpr,power=softPower) 
sizeGrWindow(10, 5)
#pdf("scale free topology power16.pdf")
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main=paste("Check Scale-free topology\n Power=",softPower))

# Calculate the adjacencies, using the soft thresholding power
softPower = 16; # Choose the lowest power from the figure #My settings
adjacency = adjacency(datExpr, power = softPower);

# Transform the adjacency into Topological Overlap Matrix (TOM), and calculate the corresponding dissimilarity
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Clustering using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()

# Set the minimum module size #My settings
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Plot the module assignment under the gene dendrogram
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
#pdf("Os_MEtree.pdf", width = 6, height = 5)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#Choose a height cut #Value of 0.25, corresponding to correlation of 0.75 
MEDissThres = 0.14 #My settings

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
#dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
# Plot the gene dendrogram again, with the original and merged module colors underneath 

#pdf("Os_Dendrogram.pdf", width = 5, height = 3)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
table(mergedColors)


# Plot the gene dendrogram after mergeing only 
#plotDendroAndColors(geneTree,mergedColors,"Merged dynamic",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)

#Save the relevant variables for use in subsequent parts if neccessary
# Rename to moduleColors
moduleColors = mergedColors


# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts if neccessary
save(MEs, moduleLabels, moduleColors, geneTree, file = "SH-orthologs.RData")

# displaying module heatmap and the eigengene barplot
person=cor(datExpr,use = 'p')
corr<-TOM
Colors<-mergedColors
colnames(corr)<-colnames(datExpr)
rownames(corr)<-colnames(datExpr)
names(Colors)<-colnames(datExpr)
colnames(person)<-colnames(datExpr)
rownames(person)<-colnames(datExpr)
umc = unique(mergedColors)
lumc = length(umc)

for (i in c(1:lumc)){
  #  if(umc[i]== "grey"){
  #    next
  #  }
  ME=MEs[, paste("ME",umc[i], sep="")]
  pdf_file_out= paste("SH-total-nostate-30-0.14-",umc[i],".pdf",sep="")
  #pdf(file = pdf_file_out, wi = 9, he = 6)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[,Colors==umc[i]])),nrgcols=30,rlabels=F,clabels=rownames(datExpr),rcols=umc[i], main=umc[i], cex.main=3)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=umc[i], main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
  #dev.off()
}


###### Step 3 ###### Exporting network data to network visualization 

# Select modules
#ReadGeneIdFiles <- grep('Sv_geneID-.*.txt', list.files(getwd()), value=TRUE)
#ReadGeneIdFiles
#FileColors <- gsub("Sv_geneID-(.*).txt", "\\1", ReadGeneIdFiles)
#modules = c(FileColors)

# Get outputs by colors
for (module in unique(mergedColors))
{
  
  # Select module probes
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, module));
  modProbes = probes[inModule];
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("SH-CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("SH-CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule]);
  
}

#Export all modules in one file
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("SH-CytoscapeInput-edges-th01.txt", sep=""),
                               nodeFile = paste("SH-CytoscapeInput-nodes-th01.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = probes,
                               altNodeNames = moduleColors);
#                               nodeAttr = moduleColors);


KIM = intramodularConnectivity(adjacency, moduleColors, scaleByMax= F)
write.csv(KIM, file = "SH-total-IntramodularConnectivity.csv")