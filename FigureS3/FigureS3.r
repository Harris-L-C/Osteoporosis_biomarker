#clusterProfiler
setwd("/Users/miaoxingrenak/Desktop/20220411文献和PPT/eg/eg4_骨质疏松/code2/OC2/code_raw data/骨质疏松数据和代码2/FigureS3/FigureS3")
options(stringsAsFactors = FALSE) 
library(WGCNA);library(clusterProfiler)
library(caret);library(reshape2)

data_g = read.csv("geneInfo.csv")
head(data_g)
dim(data_g)


colors = melt(table(data_g$moduleColor))$Var1
table(colors)
for(k in colors){
    eg = data_g[data_g$moduleColor == k,]
    enrichGO_CC <- enrichGO(gene = eg$LocusLinkID, OrgDb="org.Hs.eg.db", ont = "CC", pvalueCutoff = 1, readable= TRUE)
    enrichGO_MF <- enrichGO(gene = eg$LocusLinkID, OrgDb="org.Hs.eg.db", ont = "MF", pvalueCutoff = 1, readable= TRUE)
    enrichGO_BP <- enrichGO(gene = eg$LocusLinkID, OrgDb="org.Hs.eg.db", ont = "BP", pvalueCutoff = 1, readable= TRUE)
    enrichKEGG <- enrichKEGG(gene = eg$LocusLinkID, organism = 'hsa', pvalueCutoff = 1)


    BP_order_qvalue = enrichGO_BP@result[order(enrichGO_BP@result$qvalue),]
    BP_top_10 = BP_order_qvalue[1:10,1:9]
    BP_top_10$Type = "BP"
    BP_top_10$Module = k
    BP_top_10

    CC_order_qvalue = enrichGO_CC@result[order(enrichGO_CC@result$qvalue),]
    CC_top_10 = CC_order_qvalue[1:10,1:9]
    CC_top_10$Type = "CC"
    CC_top_10$Module = k
    CC_top_10

    MF_order_qvalue = enrichGO_MF@result[order(enrichGO_MF@result$qvalue),]
    MF_top_10 = MF_order_qvalue[1:10,1:9]
    MF_top_10$Type = "MF"
    MF_top_10$Module = k
    MF_top_10

    KEGG_order_qvalue = enrichKEGG@result[order(enrichKEGG@result$qvalue),]
    KEGG_top_10 = KEGG_order_qvalue[1:10,1:9]
    KEGG_top_10$Type = "KEGG"
    KEGG_top_10$Module = k
    KEGG_top_10

    library(treemap)
    library(ggplot2)
    library(tidyverse)
    GO_10 =rbind(rbind(rbind(BP_top_10,CC_top_10),MF_top_10),KEGG_top_10)
    a=GO_10
    GO_10$val = -log(GO_10$qvalue,10)
    godata = GO_10
    data = godata

    data$Description<-factor(data$Description,levels = unique(data$Description),ordered = T)
    Figure = ggplot(data)+geom_bar(aes(x=Description,y=-log10(pvalue),fill=Type),
    stat = 'identity')+theme(axis.text.y=element_text(vjust=1,size=15),legend.key.size = unit(0.3, "inches"),legend.title=element_text(size=15))+coord_flip()
    ggsave(Figure, file=paste(k,"function.pdf",sep="_"),width=10,height=4)
}