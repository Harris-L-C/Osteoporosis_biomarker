#Cor Module LncRNA
rm(list=ls());library(WGCNA);library(ggplot2);library(nlme)
library(gplots);library(GEOquery);library(biomaRt);library(sva);library(RColorBrewer);library(flashClust)
options(stringsAsFactors = FALSE);library(reshape2)

load("TOMcuts_2.Rdata")
traitDataT = read.csv("train_y_lncRNA.csv",row.names=1)
head(traitDataT)
datTraitsT = t(traitDataT)
head(datTraitsT)


nGenes = ncol(datExprT)
nSamples = nrow(datExprT)
datME= moduleEigengenes(datExprT, colorh1, softPower = 12)$eigengenes
datMET = orderMEs(datME)
signif(cor(datME,use="p"),2)
all.equal(rownames(datTraitsT),rownames(datME))

#heatmap
moduleTraitCorT = cor(datMET, datTraitsT, use = "p")
moduleTraitPvalueT = corPvalueStudent(moduleTraitCorT, nSamples)
moduleTraitPvalueT.FDR = p.adjust(moduleTraitPvalueT, "fdr")
dim(moduleTraitPvalueT.FDR) = dim(moduleTraitPvalueT)
dimnames(moduleTraitPvalueT.FDR) = dimnames(moduleTraitPvalueT)

pdf(file = "heatmap_moduleColorsT3.pdf",width=18,height=9)
textMatrix = paste(signif(moduleTraitCorT, 2), "\n(",signif(moduleTraitPvalueT.FDR, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorT)
par(mar = c(6, 6, 6,6))
labeledHeatmap(Matrix = moduleTraitCorT,
xLabels = colnames(datTraitsT),
yLabels = colnames(datMET),
ySymbols = colnames(datMET),cex.lab=0.60,
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.60,
zlim = c(-1,1),
main = paste("Module-LncRNA relationships"))
dev.off()