#miRNA differential analysis
rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA);library(affy);library(limma)
library(biomaRt);library(sva);library(RColorBrewer)
data.affy = ReadAffy()
data.affy
datExpr = rma(data.affy, normalize=T, background=T, verbose=T)
datExpr = exprs(datExpr)
batch = as.factor(substr(protocolData(data.affy)$ScanDate,1,10))


datMeta = read.csv("GSE63446_bene.csv",row.names=1)
head(datMeta)
datMeta$Characteristics = as.factor(datMeta$Characteristics)
datMeta$State = as.factor(datMeta$State)



RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
RNAdeg$sample.names2 = substring(RNAdeg$sample.names,1,10)
idx=  match(rownames(datMeta), RNAdeg$sample.names2)
datMeta$RNAdeg = RNAdeg$slope[idx]

datMeta$batch = as.factor(batch)
datMeta$RNAdeg = as.numeric(datMeta$RNAdeg)

datExpr[1:4,1:4]
colnames(datExpr) = substring(colnames(datExpr),1,10)
colnames(datExpr)

write.csv(datExpr,"datExpr_GSE63446.csv")
write.csv(datMeta,"datMeta_GSE63446.csv")



#Combine
library(WGCNA);library(affy);library(limma);library(biomaRt);library(sva);library(AnnotationDbi)
library(RColorBrewer);library(graph)
library(annotate);library(XML);library(IRanges);library(org.Hs.eg.db)                          
library(DBI);library(pheatmap);library(GOstats);library(base);library(edgeR)

data1 = read.csv("datExpr_GSE63446_2.csv",row.names=1,check.names=1)
data2 = read.csv("GSE64433_exprSet.csv",row.names=1,check.names=1)
data3 = read.csv("GSE91033_exprSet.csv",row.names=1,check.names=1)
dim(data1)
dim(data2)
dim(data3)
commonProbesA = intersect(intersect(rownames(data1),rownames(data2)),rownames(data3))
data1_2 = data1[commonProbesA,]
dim(data1_2)
data2_2 = data2[commonProbesA,]
dim(data2_2)
data2_3 = data3[commonProbesA,]
dim(data2_3)
save(file="shared.Rdata",data1_2,data2_2,data2_3)


all.equal(rownames(data1_2),rownames(data2_2))
all.equal(rownames(data1_2),rownames(data2_3))
dataC =cbind(cbind(data1_2,data2_2),data2_3)
dataC[1:4,1:4]
datMetaC = read.csv("dataMeta_all.csv",row.names=1,check.names=1)
dim(datMetaC)
head(datMetaC)

dataC = dataC[,rownames(datMetaC)]

idx=match(colnames(dataC), rownames(datMetaC))
datMetaC = datMetaC[idx,]
all.equal(colnames(dataC), rownames(datMetaC))

datMetaC$database = as.factor(datMetaC$database)
datMetaC$Characteristics = as.factor(datMetaC$Characteristics)
mod = model.matrix(~Characteristics, data=datMetaC)
batch1 = as.factor(datMetaC$database)
datExpr.combat = ComBat(dataC, batch=batch1, mod=mod)

plotMDS(dataC,top=20000,gene.selection="common",col=as.numeric(datMetaC$database),cex=1.8, pch=19, cex.lab = 1.8, cex.axis = 1.8)
plotMDS(datExpr.combat,top=20000,gene.selection="common",col=as.numeric(datMetaC$database),cex=1.8, pch=19, cex.lab = 1.8, cex.axis = 1.8)
plotMDS(datExpr.combat,top=20000,gene.selection="common", labels = colnames(datExpr.combat),cex.lab = 0.4,col=as.numeric(datMetaC$Characteristics),cex=1.8, pch=19, cex.lab = 1.8, cex.axis = 1.8)


#Bayes
datExpr = read.csv("datExpr_ALL_miRNA.csv",row.names=1)
datExpr[1:5,1:4]
datMeta = read.csv("dataMeta_all.csv",row.names=1)
datMeta[1:5,1:4]
#load("GSE56815_GPL96_4.RData")
all.equal(rownames(datMeta),colnames(datExpr))
Characteristics<-factor(datMeta[,"Characteristics"])       
design<-model.matrix(~-1+Characteristics)
design            
colnames(design)<-c("s1", "s2")                    
contrast.matrix<-makeContrasts(s2-s1,levels=design)   
fit <- lmFit(datExpr, design)
fit1<- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit1)
#results<-decideTests(fit2,method="global",adjust.method="BH",p.value=0.05,lfc=0.25)
#results<- topTable(fit2, coef=1, adjust.method="fdr",p.value=0.05, lfc=1, number=30000)
#summary(results)
dif<- topTable(fit2, coef=1, number=30000, adjust.method="BH", sort.by="B", resort.by="M")
write.csv(dif, "resultsdif_miRNA_3.csv")



#FigureS4A
dataset = read.csv("resultsdif_miRNA_3.csv",header = TRUE)
head(dataset)
cut_off_pvalue = 0.05
cut_off_logFC = 0.25
dataset$change = ifelse(dataset$P.Value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC > cut_off_logFC ,'Up','Down'),'Stable')
Figure = ggplot(
  dataset, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#D6604D", "#d2dae2","#4393C3"))+
  geom_vline(xintercept=c(-2.5,2.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
ggsave(Figure, file="FigureS4A.pdf")


#Figure1B
diff=read.csv("resultsdif_miRNA_3.csv",row.names=1)
head(diff)
a = read.csv("datExpr_ALL_miRNA.csv",row.names=1)
head(a)
selected = a[rownames(diff), ]
head(selected)
bene = read.csv("dataMeta_all.csv",row.names=1)
head(bene)
all.equal(colnames(selected), rownames(bene))

annotation_col = data.frame(
  Factor = as.factor(bene$Characteristics),
  RNAdeg = as.numeric(bene$RNAdeg))
head(annotation_col)
pdf(file = "FigureS4B.pdf",width = 5,height=5)
rownames(annotation_col) = colnames(selected)
pheatmap(selected, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
fontsize_row = 7, scale = "row", border_color = NA, annotation_col = annotation_col,show_rownames = F,
show_colnames = F, cluster_rows = F)
dev.off()

