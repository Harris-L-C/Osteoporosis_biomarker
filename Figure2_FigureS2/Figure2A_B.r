#GSE7429
rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA);library(affy)
library(limma);library(biomaRt)
library(sva);library(RColorBrewer);library(ggplot2);library(pheatmap)
data.affy = ReadAffy()
data.affy
datExpr = rma(data.affy, normalize=T, background=T, verbose=T)
datExpr = exprs(datExpr)
batch = as.factor(substr(protocolData(data.affy)$ScanDate,1,8))

datMeta = read.csv("GSE7429.csv",row.names=1)
head(datMeta)

RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
RNAdeg$sample.names2 = substring(RNAdeg$sample.names,1,9)
idx=  match(rownames(datMeta), RNAdeg$sample.names2)
datMeta$RNAdeg = RNAdeg$slope[idx]

datMeta$batch = as.factor(batch)
datMeta$RNAdeg = as.numeric(datMeta$RNAdeg)
head(datMeta)

#Annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")
identifier <- "affy_hg_u133a"
getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol","mirbase_id", "chromosome_name","start_position", "end_position","gene_biotype")
geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
head(geneDat) 
dim(geneDat)
idx = match(rownames(datExpr), geneDat$affy_hg_u133a)
datProbes = cbind(rownames(datExpr), geneDat[idx,])
all.equal(rownames(datMeta),colnames(datExpr))



## batch Correction
mod = model.matrix(~Characteristics, data=datMeta)
batch = as.factor(datMeta$batch)
datExpr.combat = ComBat(as.matrix(datExpr), batch=batch, mod=mod)
datExpr = datExpr.combat

# Collapse Rows
realGenes = !is.na(datProbes$external_gene_name)  #
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$external_gene_name, rowID = datProbes$affy_hg_u133a) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133a)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$external_gene_name

##Regress covariates
X = model.matrix(~Characteristics+Age+RNAdeg+batch, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
dim(X)
to_regress = (as.matrix(X[,3:15]) %*% (as.matrix(beta[3:15,])))
datExpr = datExpr - t(to_regress)

write.csv(datExpr,"datExpr_GSE7429_GPL96.csv")
write.csv(datMeta,"datMeta_GSE7429_GPL96.csv")
dim(datExpr)

#Bayes
all.equal(rownames(datMeta),colnames(datExpr))
Characteristics<-factor(datMeta[,"Characteristics"])       
design<-model.matrix(~-1+Characteristics)
design            
colnames(design)<-c("s1", "s2")                    
contrast.matrix<-makeContrasts(s1-s2,levels=design)   
fit <- lmFit(datExpr, design)
fit1<- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit1)
#results<-decideTests(fit2,method="global",adjust.method="BH",p.value=0.05,lfc=0.25)
results<- topTable(fit2, coef=1, adjust.method="fdr",p.value=0.05, lfc=1, number=30000)
summary(results)
dif<- topTable(fit2, coef=1, number=30000, adjust.method="BH", sort.by="B", resort.by="M")
dif2<-dif[dif[,"adj.P.Val"]<0.05,]
dim(dif2) 
write.csv(dif, "resultsdif_GSE7429_OC.csv")
write.csv(dif2, "resultsdif_GSE7429_OC_2.csv")



#GSE56815
rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA);library(affy)
library(limma);library(biomaRt)
library(sva);library(RColorBrewer)
data.affy = ReadAffy()
data.affy
datExpr = rma(data.affy, normalize=T, background=T, verbose=T)
datExpr = exprs(datExpr)
batch = as.factor(substr(protocolData(data.affy)$ScanDate,1,8))


datMeta = read.csv("GSE56815.csv",row.names=1)
head(datMeta)

RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
RNAdeg$sample.names2 = substring(RNAdeg$sample.names,1,10)
idx=  match(rownames(datMeta), RNAdeg$sample.names2)
datMeta$RNAdeg = RNAdeg$slope[idx]

datMeta$batch = as.factor(batch)
datMeta$RNAdeg = as.numeric(datMeta$RNAdeg)
head(datMeta)

#Annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")
identifier <- "affy_hg_u133a"
getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol","mirbase_id", "chromosome_name","start_position", "end_position","gene_biotype")
geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
head(geneDat) 
dim(geneDat)
idx = match(rownames(datExpr), geneDat$affy_hg_u133a)
datProbes = cbind(rownames(datExpr), geneDat[idx,])
colnames(datExpr) = substring(colnames(datExpr),1,10)
colnames(datExpr)
all.equal(rownames(datMeta),colnames(datExpr))

## batch Correction
mod = model.matrix(~Characteristics, data=datMeta)
batch = as.factor(datMeta$batch)
datExpr.combat = ComBat(as.matrix(datExpr), batch=batch, mod=mod)
datExpr = datExpr.combat

# Collapse Rows
realGenes = !is.na(datProbes$external_gene_name)  #
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$external_gene_name, rowID = datProbes$affy_hg_u133a) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133a)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$external_gene_name
all.equal(rownames(datMeta),colnames(datExpr))

##Regress covariates
X = model.matrix(~Characteristics+RNAdeg, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
dim(X)
to_regress = (as.matrix(X[,3:40]) %*% (as.matrix(beta[3:40,])))
datExpr = datExpr - t(to_regress)

write.csv(datExpr,"datExpr_GSE56815_GPL96.csv")
write.csv(datMeta,"datMeta_GSE56815_GPL96.csv")
dim(datExpr)

#Bayes
all.equal(rownames(datMeta),colnames(datExpr))
Characteristics<-factor(datMeta[,"Characteristics"])       
design<-model.matrix(~-1+Characteristics)
design            
colnames(design)<-c("s1", "s2")                    
contrast.matrix<-makeContrasts(s1-s2,levels=design)   
fit <- lmFit(datExpr, design)
fit1<- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit1)
#results<-decideTests(fit2,method="global",adjust.method="BH",p.value=0.05,lfc=0.25)
results<- topTable(fit2, coef=1, adjust.method="fdr",p.value=0.05, lfc=1, number=30000)
summary(results)
dif<- topTable(fit2, coef=1, number=30000, adjust.method="BH", sort.by="B", resort.by="M")
dif2<-dif[dif[,"adj.P.Val"]<0.05,]
dim(dif2) 
write.csv(dif, "resultsdif_GSE56815_OC.csv")
write.csv(dif2, "resultsdif_GSE56815_OC_2.csv")

#Combine
data1 = read.csv("datExpr_GSE56815_GPL96.csv",row.names=1,check.names=1)
data2 = read.csv("datExpr_GSE7429_GPL96.csv",row.names=1,check.names=1)
dim(data1)
dim(data2)
commonProbesA = intersect(rownames(data1),rownames(data2))
data1_2 = data1[commonProbesA,]
dim(data1_2)
data2_2 = data2[commonProbesA,]
dim(data2_2)
datExpr.combat = cbind(data1_2,data2_2)

train = data1_2
test = data2_2
train_set = read.csv("train_set_GSE56815.csv")
test_set = read.csv("test_set_GSE7429.csv")
save(file="ML_data_set.Rdata",train_set,train,test_set,test,datExpr.combat)





#Figure1A
dataset = read.csv("resultsdif_GSE7429_OC.csv",header = TRUE)
head(dataset)
cut_off_pvalue = 0.05
cut_off_logFC = 0.25
dataset$change = ifelse(dataset$adj.P.Val < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC > cut_off_logFC ,'Up','Down'),'Stable')
Figure = ggplot(
  dataset, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
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
ggsave(Figure, file="Figure1A.pdf")


#Figure1B
diff=read.csv("resultsdif_GSE56815_OC.csv",row.names=1)
head(diff)
a = read.csv("datExpr_GSE56815_GPL96.csv",row.names=1)
head(a)
selected = a[rownames(diff), ]
head(selected)
bene = read.csv("datMeta_GSE56815_GPL96.csv",row.names=1)
head(bene)
all.equal(colnames(selected), rownames(bene))

annotation_col = data.frame(
  Factor = as.factor(bene$Characteristics),
  RNAdeg = as.numeric(bene$RNAdeg))
head(annotation_col)
pdf(file = "Figure1B.pdf",width = 5,height=5)
rownames(annotation_col) = colnames(selected)
pheatmap(selected, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
fontsize_row = 7, scale = "row", border_color = NA, annotation_col = annotation_col,show_rownames = F,
show_colnames = F, cluster_rows = F, cluster_cols = F)
dev.off()