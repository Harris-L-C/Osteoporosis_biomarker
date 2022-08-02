#WGCNA
rm(list=ls());library(WGCNA);library(ggplot2);library(nlme)
library(gplots);library(GEOquery);library(biomaRt);library(sva);library(RColorBrewer);library(flashClust)
options(stringsAsFactors = FALSE);library(reshape2)
dir()

Toll = read.csv("train_data_GSE56815.csv",row.names=1)
Toll[seq(1,10),seq(1,4)]
dim(Toll)
colnames(Toll)


#芯片注释
getinfo <- c("external_gene_name","ensembl_gene_id","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="https://www.ensembl.org")
geneAnno1 <- getBM(attributes = getinfo,filters=c("external_gene_name"),values=rownames(Toll),mart=mart)
b <- geneAnno1[match(rownames(Toll),geneAnno1[,1]),]
dim(b)
head(b)
table(substr(b[,8],nchar(b[,8])-13,nchar(b[,8]))=="protein_coding")
Name<-unique(b[substr(b[,8],nchar(b[,8])-13,nchar(b[,8]))=="protein_coding",1])
Toll_1<-Toll[Name,]


datExproT = as.data.frame(t(Toll_1))
gsga = goodSamplesGenes(datExproT, verbose = 3)
gsga$allOK
gsga$goodGenes
gsga$goodSamples
if (!gsga$allOK){
datExproT = datExproT[gsga$goodSamples, gsga$goodGenes]
}
dim(datExproT)


sampleTreeT= hclust(dist(datExproT), method = "average")
pdf(file = "sampleTreeT.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeT, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

clustT = cutreeStatic(sampleTreeT, cutHeight =45, minSize = 10)
table(clustT)

keepSamples = (clustT==1)
datExprT = datExproT[keepSamples, ]
naGenes = ncol(datExprT)
naSamples = nrow(datExprT)


powers = c(seq(1, 30, by = 1))
sfta = pickSoftThreshold(datExprT,networkType="signed",corFnc="cor", powerVector = powers, verbose = 5)
pdf(file = "Figure2A_1.pdf")
par(mfrow = c(1,1))
par(mar = c(5,5,5,5))
cex1 = 0.9
plot(sfta$fitIndices[,1], -sign(sfta$fitIndices[,3])*sfta$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(sfta$fitIndices[,1], -sign(sfta$fitIndices[,3])*sfta$fitIndices[,2],
labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
abline(h=0.80,col="red")
dev.off()
pdf(file = "Figure2A_2.pdf")
par(mar = c(5,5,5,5))
plot(sfta$fitIndices[,1], sfta$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sfta$fitIndices[,1], sfta$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.0,col="red")
dev.off()




#network dendrogram
ADJ1=abs(cor(datExprT,use="p"))^12
k=softConnectivity(datE=datExprT,power=12)
datExprT=datExprT[, rank(-k,ties.method="first" )<=12044]
dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
#save(file="dissTOM.Rdata",dissTOM)

geneTree= hclust(as.dist(dissTOM),method="average")
minModSize <- 100 # Modules are at least 100 genes large
dthresh <- 0.15 # MEs are no more than 0.85 correlated, if they are then the modules are merged and the ME is re-calculated
ds <- 4 # deep split parameter to determine how finely to cut the tree
tree = cutreeDynamic(dendro = geneTree, pamStage=FALSE,minClusterSize = minModSize, cutHeight = 0.99999, deepSplit = ds, distM = as.matrix(dissTOM))


mLabelh <- "Modules"
table(tree)
dynamicColors = labels2colors(tree)
table(dynamicColors)

  
pdf(file = "net_colorh1.pdf")
plotDendroAndColors(geneTree, dynamicColors, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Defined Modules")
dev.off()


MEDissThres = dthresh
MEList = moduleEigengenes(datExprT, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");


merge = mergeCloseModules(datExprT, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;



pdf(file = "Figure2B.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

colorh1 = mergedColors
table(colorh1)
save(geneTree,colorh1,datExprT,file="TOMcuts_2.Rdata")


PLOT_PDF = T
datExpr<-datExprT
colors<-colorh1
dim(datExpr)
length(colors)

MEs=moduleEigengenes(datExpr, colors)
eigmat = MEs$eigengenes
colnames(eigmat) = gsub("ME","",colnames(eigmat))
kME = signedKME(datExpr, MEs$eigengenes)
colnames(kME) = gsub("kME", "", colnames(kME))
datMeta=read.csv("train_set_GSE56815.csv",row.names=1)
head(datMeta)
dim(datMeta)

all_colors = unique(colors)
all_colors = all_colors[!grepl("grey",all_colors)]
all_genes = colors
names(all_genes) = rownames(datExpr)


moduleTraitP = matrix(NA,nrow=length(all_colors),ncol=5)
colnames(moduleTraitP) = c("ANOVA", "Characteristics","State" ,"RNAdeg", "batch")
rownames(moduleTraitP) = all_colors
moduleTraitP

moduleTraitB = moduleTraitSE = moduleTraitP
for (m in all_colors) {
me_name = paste("ME", m, sep="")
me = MEs$eigengenes[[me_name]]
i = which(m == rownames(moduleTraitP)) 
mixedmodel = lme(me ~ Characteristics + State + RNAdeg + batch, data = datMeta, random = ~1|SetType)
moduleTraitP[i,"ANOVA"]=anova(mixedmodel)["Characteristics","p-value"] 
mixedmodel = summary(mixedmodel)$tTable
for(var in c("Characteristics","State" ,"RNAdeg")) {
moduleTraitP[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 5]
moduleTraitB[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 1]
moduleTraitSE[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 2]
}
for(cov in c("Characteristics","State" ,"RNAdeg", "batch")) {
mod = summary(lm(me ~ datMeta[,cov]))
moduleTraitP[i,cov] = mod$coefficients[2, 4]
moduleTraitB[i,cov] = mod$coefficients[2, 1]
moduleTraitSE[i,cov] = mod$coefficients[2, 2]
}
  }

moduleTraitP.fdr = p.adjust(moduleTraitP, "fdr")
dim(moduleTraitP.fdr) = dim(moduleTraitP)
dimnames(moduleTraitP.fdr) = dimnames(moduleTraitP)
out = cbind(moduleTraitB, moduleTraitP, moduleTraitP.fdr)
colnames(out)[1:5] = paste("Beta.", colnames(out)[1:5], sep="")
colnames(out)[6:10] = paste("P.", colnames(out)[6:10], sep="")
colnames(out)[11:15] = paste("FDR.", colnames(out)[11:15], sep="")

row_idx = rownames(moduleTraitB)
col_idx = c(2:5)
bpdata = melt(moduleTraitB[row_idx,col_idx])
semdata = melt(moduleTraitSE[row_idx,col_idx])
pdata = melt(moduleTraitP.fdr[row_idx,col_idx])
bpdata$sem = semdata$value
bpdata$p = pdata$value
bpdata$p.symbol = ""
bpdata$p.symbol[bpdata$p<0.05] = "*"
bpdata$p.symbol[bpdata$p<0.01] = "**"
bpdata$p.symbol[bpdata$p<0.001] = "***"

bpdata$X1 = factor(bpdata$X1, levels=rownames(moduleTraitB))
bpdata
dim(bpdata)

bpdata2<-bpdata[c(1:10),]
Fig3C.byCharacteristics2= ggplot(bpdata2, aes(x=Var1, y=value,fill=Var1,group=Var1, label=p.symbol))+ facet_wrap(~Var2,ncol=1) + 
geom_bar(stat="identity", position=position_dodge(), color="black") +
geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.35,width=.35) +
theme_minimal() + scale_fill_manual(name="Factor",values=levels(bpdata2$Var1)) +
labs(y="BETA", x="") +
geom_text(color="red",size=2.5,aes(y=value+ sign(value)*sem + sign(value)*.01, angle = 90), position=position_dodge(.9))  + 
scale_x_discrete() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"),legend.position = "none",
axis.text.y = element_text(size=8),
axis.text.x = element_text(angle=45, hjust=1))+coord_flip()

pdf(file = "Figure2C.pdf",width = 2.5, height = 3)
par(mar = c(5,5,8,8))
Fig3C.byCharacteristics2
dev.off()



#load("TOMcuts_2.Rdata")
traitDataT = read.csv("train_set_GSE56815.csv",row.names=1)
head(traitDataT)
all.equal(rownames(traitDataT),rownames(datExprT))
nGenes = ncol(datExprT)
nSamples = nrow(datExprT)
datME= moduleEigengenes(datExprT, colorh1, softPower = 12)$eigengenes
datMET = orderMEs(datME)
signif(cor(datME,use="p"),2)
Factor = as.data.frame(traitDataT$Factors);
names(Factor) = "Factor"
modNames = substring(names(datMET), 3)

geneModuleMembership = as.data.frame(cor(datExprT, datMET, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExprT, Factor, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Factor), sep="");
names(GSPvalue) = paste("p.GS.", names(Factor), sep="");



annot = read.csv(file = "datProbes_GPL96.csv")
dim(annot)
names(annot)
probes = names(datExprT)                                  #colnames or names
probes2annot = match(probes, annot$external_gene_name)
sum(is.na(probes2annot))
geneInfo0 = data.frame(substanceBXH = probes,
                    geneSymbol = annot$external_gene_name[probes2annot],
                    LocusLinkID = annot$entrezgene_id[probes2annot],
                    moduleColor = colorh1,
                    geneTraitSignificance,
                    GSPvalue)

modOrder = order(-abs(cor(datMET, Factor, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                       MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                     paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Factor));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")