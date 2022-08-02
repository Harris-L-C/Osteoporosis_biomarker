#LASSO
#load("ML_data_set.Rdata")
library(ggplot2)
library(ROCR)
library(glmnet)
library(pROC)
library(gplots)
library(reshape2)

module = read.csv("top25.csv")
head(module)


train_y <- train[module$X,]
test_y <-  test[module$X,]

train_x <- train_set
test_x <-  test_set

CoxExample = list(x=t(train_y),y=as.numeric(as.factor(train_x$Characteristics)))

x <- CoxExample$x
y <- CoxExample$y

pdf(file = "Figure1C_1.pdf",width = 5,height=5)
fit <- glmnet(as.matrix(x), y,family ="gaussian")
plot(fit)
dev.off()

pdf(file = "Figure1C_2.pdf",width = 5,height=5)
cvfit <- cv.glmnet(as.matrix(x), y, family = "gaussian", type.measure = "auc")
plot(cvfit)
dev.off()


library(pROC)
library(gplots)

preds <- predict(cvfit, newx = as.matrix(t(test_y)), type = 'response')
perf <- performance(prediction(preds, as.numeric(as.factor(test_x$Characteristics))), 'tpr', 'fpr')
roc_curve <- roc(as.numeric(as.factor(test_x$Characteristics)),preds)

a = plot(roc_curve, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightgreen", print.thres=TRUE,main="Test Data ROC curve")
ggsave(a, file="Figure1F.pdf") 

preds_tra <- predict(cvfit, newx = as.matrix(t(train_y)), type = 'response')
perf_tra <- performance(prediction(preds_tra, as.numeric(as.factor(train_x$Characteristics))), 'tpr', 'fpr')
roc_curve_tra<- roc(as.numeric(as.factor(train_x$Characteristics)),preds_tra)

a = plot(roc_curve_tra, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="red", print.thres=TRUE,main="Train Data ROC curve")
ggsave(a, file="Figure1E.pdf") 


#Figure1D
iris_input_test = test_y
head(iris_input_test)
iris_input_datMeta = test_x
head(iris_input_datMeta)
all.equal(rownames(iris_input_datMeta),colnames(iris_input_test))
pca1 <- prcomp(t(iris_input_test),center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = iris_input_datMeta$Characteristics))+
  stat_ellipse(aes(fill = iris_input_datMeta$Characteristics),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("lightblue","tan"))+
  scale_colour_manual(values = c("lightblue","tan"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(p.pca1, file="Figure1D.pdf") 



#FigureS2
for(i in c(1:(length(colnames(t(train_y)))))){
   rocobj<-roc(as.numeric(as.factor(train_x$Characteristics)),t(train_y)[,i])
   print(rocobj)
   pdf(file = paste(colnames(t(train_y))[i],"ROC.pdf",sep="_"),width = 5,height=5)
   a = plot(rocobj, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
   grid.col=c("blue", "pink"), max.auc.polygon=TRUE,
   auc.polygon.col="red", print.thres=TRUE,main=paste(colnames(t(train_y))[i],"ROC curve",sep=" "))
   dev.off()   
}


for(j in c(1:(length(colnames(t(test_y)))))){
   rocobj<-roc(as.numeric(as.factor(test_x$Characteristics)),t(test_y)[,j])
   print(rocobj)
   pdf(file = paste(colnames(t(test_y))[j],"ROC.pdf",sep="_"),width = 5,height=5)
   a = plot(rocobj, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
   grid.col=c("blue", "pink"), max.auc.polygon=TRUE,
   auc.polygon.col="lightgreen", print.thres=TRUE,main=paste(colnames(t(test_y))[j],"ROC curve",sep=" "))
   dev.off()   
}