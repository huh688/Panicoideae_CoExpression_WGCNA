setwd("~/Downloads/WGCNA_YQ/Module similarity")

library("WGCNA")

load("/Users/hao/Downloads/WGCNA/Bd/Bd_min30cut14s16/Bd-orthologs.RData")
BdLabels = moduleLabels;
BdColors = moduleColors;
BdTree = geneTree;
BdMEs = orderMEs(MEs, greyName = "ME0");

load("/Users/hao/Downloads/WGCNA/Os/Os_min30cut14s16/SH-orthologs.RData")

SHLabels = moduleLabels;
SHColors = moduleColors;
SHTree = geneTree;
SHMEs = orderMEs(MEs, greyName = "ME0");

# Isolate the module labels in the order they appear in ordered module eigengenes
BdModuleLabels = substring(names(BdMEs), 3)
SHModuleLabels = substring(names(SHMEs), 3)
# Convert the numeric module labels to color labels
BdModules = unique(BdColors)
SHModules = unique(SHColors)
# Numbers of female and consensus modules
nBdMods = length(BdModules)
nSHMods = length(SHModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nBdMods, ncol = nSHMods);
CountTbl = matrix(0, nrow = nBdMods, ncol = nSHMods);
# Execute all pairwaise comparisons - to produce p-value from fisher's test
for (fmod in 1:nBdMods)
  for (cmod in 1:nSHMods)
  {
    BdMembers = (BdColors == BdModules[fmod]);
    SHMembers = (SHColors == SHModules[cmod]);
    pTable[fmod, cmod] = fisher.test(BdMembers, SHMembers, alternative = "greater")$p.value;
    CountTbl[fmod, cmod] = sum(BdColors == BdModules[fmod] & SHColors ==
                                 SHModules[cmod])
  }

rownames(pTable) = paste("Bd", BdModules)
colnames(pTable) = paste("SH", SHModules)
write.csv(pTable, file = "pTable-BdtoSH.csv")

# Execute all pairwaise comparisons again - to produce -log10 transfromed p-value
for (fmod in 1:nBdMods)
  for (cmod in 1:nSHMods)
  {
    BdMembers = (BdColors == BdModules[fmod]);
    SHMembers = (SHColors == SHModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(BdMembers, SHMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(BdColors == BdModules[fmod] & SHColors ==
                                 SHModules[cmod])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
BdModTotals = apply(CountTbl, 1, sum)
SHModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(9, 13, 4, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", SHModules),
               yLabels = paste(" ", BdModules),
               colorLabels = TRUE,
               xSymbols = paste("SH ", SHModules, ": ", SHModTotals, sep=""),
               ySymbols = paste("Bd ", BdModules, ": ", BdModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Module similarity between Bd and SH orthologs subnetworks",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
#dev.off();

write.csv(pTable, file = "pTable-log-BdtoSH.csv")

