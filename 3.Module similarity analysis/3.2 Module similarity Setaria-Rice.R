setwd("~/Downloads/WGCNA_YQ/Module similarity")

library("WGCNA")

load("/Users/hao/Downloads/WGCNA/Sv/Sv_min30cut14s16/ME34-orthologs.RData")
ME34Labels = moduleLabels;
ME34Colors = moduleColors;
ME34Tree = geneTree;
ME34MEs = orderMEs(MEs, greyName = "ME0");

load("/Users/hao/Downloads/WGCNA/Os/Os_min30cut14s16/SH-orthologs.RData")

SHLabels = moduleLabels;
SHColors = moduleColors;
SHTree = geneTree;
SHMEs = orderMEs(MEs, greyName = "ME0");

# Isolate the module labels in the order they appear in ordered module eigengenes
ME34ModuleLabels = substring(names(ME34MEs), 3)
SHModuleLabels = substring(names(SHMEs), 3)
# Convert the numeric module labels to color labels
ME34Modules = unique(ME34Colors)
SHModules = unique(SHColors)
# Numbers of female and consensus modules
nME34Mods = length(ME34Modules)
nSHMods = length(SHModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nME34Mods, ncol = nSHMods);
CountTbl = matrix(0, nrow = nME34Mods, ncol = nSHMods);
# Execute all pairwaise comparisons - to produce p-value from fisher's test
for (fmod in 1:nME34Mods)
  for (cmod in 1:nSHMods)
  {
    ME34Members = (ME34Colors == ME34Modules[fmod]);
    SHMembers = (SHColors == SHModules[cmod]);
    pTable[fmod, cmod] = fisher.test(ME34Members, SHMembers, alternative = "greater")$p.value;
    CountTbl[fmod, cmod] = sum(ME34Colors == ME34Modules[fmod] & SHColors ==
                                 SHModules[cmod])
  }

rownames(pTable) = paste("ME34", ME34Modules)
colnames(pTable) = paste("SH", SHModules)
write.csv(pTable, file = "pTable-ME34toSH.csv")

# Execute all pairwaise comparisons again - to produce -log10 transfromed p-value
for (fmod in 1:nME34Mods)
  for (cmod in 1:nSHMods)
  {
    ME34Members = (ME34Colors == ME34Modules[fmod]);
    SHMembers = (SHColors == SHModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(ME34Members, SHMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(ME34Colors == ME34Modules[fmod] & SHColors ==
                                 SHModules[cmod])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
ME34ModTotals = apply(CountTbl, 1, sum)
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
               yLabels = paste(" ", ME34Modules),
               colorLabels = TRUE,
               xSymbols = paste("SH ", SHModules, ": ", SHModTotals, sep=""),
               ySymbols = paste("ME34 ", ME34Modules, ": ", ME34ModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Module similarity between ME34 and SH orthologs subnetworks",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
#dev.off();

write.csv(pTable, file = "pTable-log-ME34toSH.csv")

