






#source-packages
library(tidyverse)
library(glmnet)
library(ggplot2)
library(survival)
#########################################
#'load-filtered-features
(load("1-OS-and-filtered-signatures-matrix.RData"))
############################################
sig_eset[1:7,1:8]
##################################
##################################
#bootstraping combine with LASSO cox regression
##################################
#10000 iteration
##################################
cao<-as.character()
for(i in 1:10000){
  a<-floor(runif(floor(dim(sig_eset)[1])*0.80,1,dim(sig_eset)[1]))
  rt<-as.matrix(sig_eset[a,1:2])
  genedatasd<-as.matrix(as.data.frame(sig_eset[a,3:ncol(sig_eset)]))
  fit2<-cv.glmnet(genedatasd, rt, family="cox", maxit=1000,nfolds = 10,alpha=1)
  myCoefs <- coef(fit2, s="lambda.min");
  b<-myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
  cao<-append(cao,b)
}
################################
aa<-as.data.frame(sort(table(cao),decreasing = T))
################################


