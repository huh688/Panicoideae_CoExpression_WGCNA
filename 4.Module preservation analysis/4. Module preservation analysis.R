library("WGCNA")

options(stringsAsFactors = FALSE);
allowWGCNAThreads()

setwd("~/Downloads/WGCNA/concensus")

DataSet = read.csv("/Users/hao/Downloads/WGCNA_YQ/orthologs_expressed_counts.csv", header=T)


# We work with two sets:
nSets = 3;

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Sv", "Bd", "Os")

# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = list()

multiExpr[[1]] = list(data = as.data.frame(t(DataSet[, -c(1:41)])));
#names(multiExpr[[1]]$data) = DataSet$Sv_geneID;
#rownames(multiExpr[[1]]$data) = names(DataSet[, -c(1:41)]);

multiExpr[[2]] = list(data = as.data.frame(t(DataSet[, -c(1:5, 24:59)])));
#names(multiExpr[[2]]$data) = DataSet$Bd_geneID;
#rownames(multiExpr[[2]]$data) = names(DataSet[, -c(1:5, 24:59)]);


multiExpr[[3]] = list(data = as.data.frame(t(DataSet[, -c(1:23, 42:59)])));
#names(multiExpr[[3]]$data) = DataSet$Os_geneID;
#rownames(multiExpr[[3]]$data) = names(DataSet[, -c(1:23, 42:59)]);

names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)

# Loading of module labels
load("/Users/hao/Downloads/WGCNA/Sv/Sv_min30cut14s16/ME34-orthologs.RData")
SvLabels = moduleLabels;
SvColors = moduleColors;
SvTree = geneTree;
SvMEs = orderMEs(MEs, greyName = "ME0");

load("/Users/hao/Downloads/WGCNA/Bd/Bd_min30cut14s16/Bd-orthologs.RData")
BdLabels = moduleLabels;
BdColors = moduleColors;
BdTree = geneTree;
BdMEs = orderMEs(MEs, greyName = "ME0");

load("/Users/hao/Downloads/WGCNA/Os/Os_min30cut14s16/SH-orthologs.RData")
OsLabels = moduleLabels;
OsColors = moduleColors;
OsTree = geneTree;
OsMEs = orderMEs(MEs, greyName = "ME0");

# Create an object (list) holding the module labels for each set:
colorList = list(SvColors, BdColors, OsColors);
 
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels;


load("/Users/hao/Downloads/WGCNA_YQ/concensus/SvBdOs-modulePreservation-max3000.RData")

ref = 3 # Select the human data as reference
test = 2 # Select the chimp data as test
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

stats= mp$preservation$observed[[ref]][[test]];
labelsX = rownames(stats)
labelsX[labelsX=="gold"] = "gold"
modColors = labelsX;
plotMods = !(modColors %in% c("grey", "gold"));
moduleSizes = stats[plotMods, 1];
textLabels = match(modColors, standardColors())[plotMods];
colorLabels = labelsX[plotMods];


#  prepare a list containing the x and y coordinates of scatterplots and set up necessary options.
nModules = sum(plotMods);
nPlots = 6
plotData = list();
# Fill up the plotData
plotData[[1]] = plotData[[2]] = matrix(0, nModules, nPlots);
plotData[[1]][, c(1:4)] = moduleSizes;
plotData[[2]][, 1] = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods];
plotData[[2]][, 2] = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods];
# Match the modulePreservation ordering of modules to that of clusterRepro
crLabels = sort(unique(colorLabels));
mp2cr = match(colorLabels, crLabels);
# Scatterplots of IGP and p-value vs. module size
plotData[[2]][, 3] = cr[[ref]][[test]]$Actual.IGP
plotData[[2]][, 4] = -log10(cr[[ref]][[test]]$p.value + 1e-4);
# Scatterplot of observed IGP vs. Zsummary and medianRank
plotData[[1]][, c(5,6)] = plotData[[2]][, c(1:2)];
plotData[[2]][, c(5,6)] = plotData[[2]][, 3];
# Plot annotation
xLabs = c(rep("Module size", 4), "Zsummary", "Median rank");
yLabs = c("Zsummary", "Median rank", "Observed IGP", "-log10(IGP perm p)", "Observed IGP", "Observed IGP");
mains = spaste(LETTERS[1:nPlots], ". ", #rep("Ref: Human, Test: Chimp\n", nPlots),
               c(yLabs[1:4], paste(yLabs[5:6], "vs.", xLabs[5:6])),
               c("", "", "", "", "\n", "\n"));
# Scatterplot options
verbose = c(rep(FALSE, 4), rep(TRUE, 2));
ablines = list(c(0, 2, 10), NA, NA, c(-log10(0.05), -log10(0.05/nModules)), NA, NA);
abColors = list(c("black", "blue", "darkgreen"), NA, NA, c("blue", "red"), NA, NA);
logs = c("x", "x", "x", "x", "", "");
invertY = c(FALSE, TRUE, rep(FALSE, 4));
verSP = function(...) { verboseScatterplot(..., abline = TRUE) }


# The actual plotting code starts here.
cexLabels = 1.4
sizeGrWindow(6,9);
#pdf(file = "Plots/HumanSpecific-NetworkAndIGPStatistics.pdf", w=6, h=9, onefile = FALSE);
par(mfrow = c(3,2));
par(mar = c(3.3, 3.3, 3.2, 0.5));
par(mgp = c(2, 0.6, 0))
for (p in 1:2)
{
  x = plotData[[1]][, p];
  y = plotData[[2]][, p]
  miny = min(y, ablines[[p]], na.rm = TRUE);
  maxy = max(y, ablines[[p]], na.rm = TRUE);
  miny = miny - (maxy-miny)*0.1;
  maxy = maxy + (maxy-miny)*0.1;
  (if (verbose[p]) verSP else plot ) (plotData[[1]][, p], plotData[[2]][, p],
                                      main = mains[p],
                                      xlab = xLabs[p],
                                      ylab = yLabs[p],
                                      cex.main = cexLabels, cex.lab = cexLabels, cex.axis = cexLabels,
                                      bg = colorLabels,
                                      col = colorLabels, cex = 2.2,
                                      ylim = if (invertY[p]) c(maxy, miny) else c(miny, maxy),
                                      pch = 21,
                                      log = logs[p]);
  #labelPoints(plotData[[1]][, p], plotData[[2]][, p], textLabels, cex = cexLabels, offs = 0.02);
  labelPoints(plotData[[1]][, p], plotData[[2]][, p], colorLabels, cex = 1, offs = 0.1);
  if (!is.na(ablines[[p]][[1]]))
    for (al in 1:length(ablines[[p]]))
      abline(h = ablines[[p]][[al]], col = abColors[[p]][[al]], lty = 2);
}
# If plotting into a pdf file, close the file. An un-closed pdf file is not readable.
#dev.off();

