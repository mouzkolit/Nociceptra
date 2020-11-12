rm(list=ls())
library(WGCNA)
library("biomaRt")
library(dplyr)
options(stringsAsFactors = FALSE);
directory = "/home/physiologie/Desktop/iPSC/mRNA/RawData/analysis/counts/Network_construction"
setwd(directory)

ad3 = read.csv("AD3_miRNAs_timeline_ltr-results-with-normalized.csv") 
ad2 = read.csv("AD2_miRNAs_timeline_ltr-results-with-normalized.csv")
eight = read.csv("840_miRNAs_timeline_ltr-results-with-normalized.csv")


nSets = 3;
setLabels = c("AD3","AD2","840")
shortLabels = c("AD3","AD2","840")
multiExpr = vector(mode = "list", length = 3)
multiExpr[[1]] = list(data = as.data.frame(t(ad3[-c(1:8)])));
names(multiExpr[[1]]$data) = ad3$Row.names;
rownames(multiExpr[[1]]$data) = names(ad3)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(ad2[-c(1:8)])));
names(multiExpr[[2]]$data) = ad2$Row.names;
rownames(multiExpr[[2]]$data) = names(ad2)[-c(1:8)];
multiExpr[[3]] = list(data = as.data.frame(t(eight[-c(1:8)])));
names(multiExpr[[3]]$data) = eight$Row.namess;
rownames(multiExpr[[3]]$data) = names(eight)[-c(1:8)];


exprSize = checkSets(multiExpr)
gsg = goodSamplesGenesMS(multiExpr, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:3)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "mirna_SampleClustering_filtering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:3)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

# Network construction--> first consensus

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = 3);
# Call the network topology analysis function for each set in turn
for (set in 1:3)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]],networkType = "signed");
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:3)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:3)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}


net = blockwiseConsensusModules(
  multiExpr, power = 7, maxBlockSize = 37000, networkType = "signed",  minModuleSize = 10, deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.01, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 1)


consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];

sizeGrWindow(8,6);
pdf(file = "mirna_ConsensusDendrogram-auto_filtered.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

save(consMEs, moduleLabels, moduleColors, consTree, file = "mirna_Consensus-NetworkConstruction-auto_filterd.RData")

consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)

sizeGrWindow(8,10);
pdf(file = "mirna_EigengeneNetworks_filtered.pdf", width= 40, height = 50);
par(cex = 0.9)
plotEigengeneNetworks(consMEsC, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off();


# try to get the genes exclude the grey modul
probes = names(multiExpr[[1]]$data)
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
kME = list();
for (set in 1:nSets)
{
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z+kME[[3]]$Z)/sqrt(2);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor,kME[[3]]$cor, kME[[1]]$p, kME[[2]]$p,kME[[3]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(1647, 8*nMEs)


colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.","kME.set3", "p.kME.set1.", "p.kME.set2.","p.kME.set3", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(8, nMEs)))
info = data.frame(Probe = probes,
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  kMEmat);

write.csv(info, file = "mirna_consensusAnalysis-CombinedNetworkResults_filtered.csv",
          row.names = FALSE, quote = FALSE);
save(consMEs, moduleLabels, moduleColors, consTree, file = "mirna_End_Consensus-NetworkConstruction-auto_filtered.RData")
