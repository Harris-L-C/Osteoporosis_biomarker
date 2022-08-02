library(ggplot2)
library(ROCR)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
module = read.csv("check_genes2.csv")
head(module)
dat = read.csv("datExpr_GSE12267.csv",row.names=1)
head(dat)
dat2 <- dat[module$X,]
head(dat2)
dat3 = read.csv("datMeta_GSE12267_2.csv",row.names=1)
head(dat3)

#Figure7A_E
idx =match(rownames(dat3),colnames(dat2))
datExpr2=dat2[,idx]
all.equal(rownames(dat3),colnames(datExpr2))
datExpr=datExpr2
dat4 = cbind(t(datExpr),dat3)
colnames(dat4)
for(k in colnames(dat4[,c(1:8)])){
    print(k)
    p <- ggboxplot(dat4, x = "Stage", y = k,
               color = "Stage", palette = "aaas",
               add = "jitter")
    p2 = p + stat_compare_means(label = "p.format")
    ggsave(p2, file=paste(k,"Stage.pdf",sep="_"),width = 4.5,height=4.5) 
}


#Figure7F
dat= read.csv("datExpr_GSE12267.csv",row.names=1)
head(dat)
lncRNA_genes = read.csv("check_genes3.csv",row.names=1,check.names=T)
lncRNA = dat[rownames(lncRNA_genes),]
head(lncRNA)
lncRNA= t(lncRNA)
head(lncRNA)
Figure_drug = ggscatter(as.data.frame(lncRNA), x = "LINC00623", y = "TNFRSF11B",add = "reg.line", 
  conf.int = TRUE,add.params = list(fill = "lightblue"))+ 
  stat_cor(method = "spearman", label.x = 7.0, label.y = c(10,7.0))
ggsave(Figure_drug, file="TNFRSF11B.pdf",width = 3.5,height = 3.5) 