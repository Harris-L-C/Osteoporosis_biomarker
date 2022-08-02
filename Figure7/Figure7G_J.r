#CEL_GPL570
rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(RColorBrewer)
library(preprocessCore)
datExpr = read.csv("datExpr_GSE176086.csv",row.names=1)
head(datExpr)
genes = read.csv("genes2.csv",row.names=1)
head(genes)
datMeta = read.csv("GSE176086_pdata.csv",row.names=1)
head(datMeta)
dat = t(datExpr)
dat[1:4,1:4]
dat2 = dat[,rownames(genes)]
head(dat2)
all.equal(rownames(datMeta),rownames(dat2))
library(reshape2)
library(ggpubr)
for(k in rownames(genes)){
    print(k)
    data <- data.frame(datMeta$group2, dat2[,k])
    colnames(data) = c("group",k)
    datag = melt(data)
    datag$group <- factor(datag$group,levels = c('control','TNFa (100ng/mL)'))
    Figure = ggplot(datag, aes(x=group, y=value,color=group)) + theme(panel.background = element_blank(),axis.line = element_line())+
    geom_boxplot()+ scale_color_manual(values=c("blue", "red"))+stat_compare_means(
    comparisons = list(c('control','TNFa (100ng/mL)')),method = "wilcox.test")+labs(title = k)             
    ggsave(Figure, file=paste(k,"group2.pdf",sep="_"),width = 6.35,height=5) 
}
