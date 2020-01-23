






#source-packages
library(tidyverse)
library(glmnet)
library(ggplot2)
library(survival)
#########################################
#'load-filtered-features
(load(paste0(path_training,"1-features.RData")))
###########################################
(load(paste0(path_training,"2-pdata-target-signature-matrix.RData")))
############################################
target_eset[1:7,1:8]
target_eset<-column_to_rownames(target_eset,var = "ID")
######################################
esetlasso<-target_eset[,c("status","time",fea_names)]
######################################
######################################
dim(esetlasso);colnames(esetlasso)[1:50]
##################################
#bootstraping combine with LASSO cox regression
##################################
#10000 iteration
##################################
cao<-as.character()
for(i in 1:10000){
  a<-floor(runif(floor(dim(esetlasso)[1])*0.80,1,dim(esetlasso)[1]))
  rt<-as.matrix(esetlasso[a,1:2])
  genedatasd<-as.matrix(as.data.frame(esetlasso[a,3:ncol(esetlasso)]))
  fit2<-cv.glmnet(genedatasd, rt, family="cox", maxit=1000,nfolds = 10,alpha=1)
  myCoefs <- coef(fit2, s="lambda.min");
  b<-myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
  cao<-append(cao,b)
}
################################
aa<-as.data.frame(sort(table(cao),decreasing = T))
################################


