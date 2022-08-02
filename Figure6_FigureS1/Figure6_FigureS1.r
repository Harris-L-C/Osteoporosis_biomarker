library(pheatmap)
library(ggplot2)
library(ggpubr)

diff=read.csv("lncRNA.csv",row.names=1)
head(diff)
a<-read.csv("datExpr_GSE56815_GPL96.csv",row.names=1)
head(a)
selected <- a[rownames(diff),]
head(selected)
bene<-read.csv("datMeta_GSE56815_GPL96.csv",row.names=1)
head(bene)
all.equal(colnames(selected), rownames(bene))

#Figure6B_E
selected_all = cbind(t(selected),bene)
selected_all[1:4,1:4]
colnames(selected_all)
data = selected_all[,c(1:5)]
head(data)
colnames(data)
selected_all$Characteristics = as.factor(selected_all$Characteristics)
pdf("boxplot_State_mRNA.pdf",height=6,width=8)           
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(data))
y=c(1:ncol(data))
plot(x,y,
     xlim=c(-1,15),ylim=c(min(data),max(data)+1),cex.axis=1.2,
     main="",xlab="", ylab="Fraction",cex.lab=1,
     pch=21,
     col="white",
     xaxt="n")
legend("topright", c("low_BMD","high_BMD" ), col=c("lightblue","pink"), pch=15, cex=0.7)
for(i in 1:ncol(data)){
  normalData=data[selected_all$Characteristics == "low_BMD",i]
  tumorData=data[selected_all$Characteristics == "high_BMD",i]
  boxplot(normalData,at=3*(i-1),lty=1,add = T,col = 'lightblue', axes=FALSE,width=5,outline=F)
  boxplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'pink', axes=FALSE,width=5,outline=F)
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  #lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.1,labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
  text(seq(1,15,3),3.85,xpd = NA,labels=colnames(data),cex = 0.8,srt = 45,pos=2)
}
dev.off()




#FigureS1
library(ROCR);library(readxl);
library(caret);library(glmnet);
library(corrplot)library(Metrics);
library(ggplot2);library(pROC);
library(survival);library(gplots)
library(survminer);library(rms)
library(reshape2)
lncRNA2 = read.csv("lncRNA2_cph.csv",row.names=1)
lncRNA2$risk<-ifelse(lncRNA2$riskscore>median(lncRNA2$riskscore),"high","low")
dim(lncRNA2)
colnames(lncRNA2)
head(lncRNA2)

dir.create("risk")
setwd("risk")
for(j in colnames(lncRNA2[,c(1:25)])){
    print(j)
    p <- ggboxplot(lncRNA2, x = "risk", y = j,
               color = "risk", palette = "aaas",
               add = "jitter")
    p2 = p + stat_compare_means(label = "p.format")
    ggsave(p2, file=paste(j,"risk_level.pdf",sep="_"),width = 3.5,height=3.5) 
}