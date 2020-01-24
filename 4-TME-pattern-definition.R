
#loading packages
library(ConsensusClusterPlus)
library(pheatmap)
library(survival)
library(survminer)
library(tidyverse)
#################################
#loading-phenotype-data
(load(paste0("1-pdata-IMvigor210.RData")))
pdata<-rownames_to_column(pdata,var = "ID")
##################################
#load-signature-matrix-data
(load(paste0(path_tme,"2-TME-Cell-fraction-7-log2TPM-symbol-IMvigor210.Rdata")))
tme_combine[1:5,1:5]
#################################
#extracting-the-result-of-CIBERSORT
cibersort<-tme_combine[,c("ID",colnames(tme_combine)[grepl(colnames(tme_combine),pattern = "_cibersort")])]
cibersort<-column_to_rownames(cibersort,var = "ID")
cibersort<-scale(cibersort,center = T,scale = T)
##################################
#cunsensu-cluster
##################################
title=tempdir();tempdir()
set.seed(823924352)
#########注意这部分-需要把矩阵文件转过来
results = ConsensusClusterPlus(t(cibersort),maxK=6,reps=1000,
                               pItem=0.8,pFeature=1,title=title,
                               clusterAlg="kmdist",distance="spearman",
                               seed=1262118388.71279,plot="png")
##################################
#classifying patients using kmeans algorithm
set.seed(122437)
##################################
km<- kmeans(cibersort,2) # determin how many cluster you want, I specify 2 here
kmeans<- cbind(t(cibersort), km$cluster) # combine the cluster with the matrix
o<- order(kmeans[,ncol(kmeans)]) # order the last column
kmeans<- kmeans[order(kmeans[,ncol(kmeans)]),] # order the matrix according to the order of the last column
colnames(kmeans)[ncol(kmeans)]<-"groups"
summary(as.factor(kmeans[,ncol(kmeans)]))  
###################################
###################################
groups<-kmeans[,ncol(kmeans)]
annotation<-data.frame(TMEcluster=groups)
annotation$TMEcluster<-as.factor(as.character(annotation$TMEcluster))
annotation$TMEcluster<-revalue(annotation$TMEcluster,c("1" = "TMEclusterA","2" = "TMEclusterB"))
head(annotation)
###################################
colnames(pdata)
clin<-pdata[pdata$ID%in%rownames(annotation),]
####################################
clin$binaryResponse<-as.character(clin$binaryResponse)
clin$binaryResponse[is.na(clin$binaryResponse)]<-"NE"
clin$binaryResponse<-revalue(clin$binaryResponse,c("CR/PR"="CRPR",
                                                   "SD/PD"="SDPD",
                                                   "NE"="NE"))
##################################
clin$IC_Level<-as.character(clin$IC_Level)
clin$IC_Level[is.na(clin$IC_Level)]<-"NE"
clin$IC_Level<-revalue(clin$IC_Level,c("IC0"="IC0","IC1"="IC1","IC2+"="IC2","NE"="NE"))
summary(as.factor(clin$IC_Level))
####################################
clin$TC_Level<-as.character(clin$TC_Level)
clin$TC_Level[is.na(clin$TC_Level)]<-"NE"
clin$TC_Level<-revalue(clin$TC_Level,c("TC0"="TC0","TC1"="TC1","TC2+"="TC2","NE"="NE"))
summary(as.factor(clin$TC_Level))
#################################
clin$Immune_phenotype<-as.character(clin$Immune_phenotype)
clin$Immune_phenotype[is.na(clin$Immune_phenotype)]<-"NE"
summary(as.factor(clin$Immune_phenotype))
####################################
clin<-clin[,c("ID","Best_Confirmed_Overall_Response","binaryResponse","IC_Level","TC_Level",
              "Immune_phenotype","Tobacco_Use_History")]
########################################
annotation<-rownames_to_column(annotation,var = "ID")
annotation<-merge(annotation,clin,by="ID")
head(annotation)
colnames(annotation)[3]<-"BOR"
colnames(annotation)[8]<-"Tobacco_use"
annotation<-column_to_rownames(annotation,var = "ID")
#############################################
sum(is.na(annotation))
colnames(annotation)
###############################################
ann_colors = list(
  TMEcluster = c(TMEclusterA = "#FDBF6F",TMEclusterB= "#1F78B4"),
  BOR=c(NE="#666666",CR="#1F78B4",PR="#33A02C",SD="#E6AB02",PD="#CAB2D6"),
  binaryResponse=c(NE= "#666666",CRPR="#B2DF8A",SDPD="#CAB2D6"),
  IC_Level=c(IC0="#D1E5F0",IC1="#92C5DE",IC2="#4393C3",NE="#666666"),
  TC_Level=c(TC0="#D1E5F0",TC1="#92C5DE",TC2="#4393C3",NE="#666666"),
  Immune_phenotype=c(excluded="#92C5DE",desert="#FDBF6F",inflamed="#33A02C",NE="#666666"),
  Tobacco_use=c(PREVIOUS= "#FB9A99",NEVER="#1B9E77",CURRENT="#67001F"))  ###  "#FF7F00" "#1F78B4","#FDBF6F","#FF7F00","#D1E5F0"
colnames(kmeans)
###############################################
##############################################
for (i in 1:c(ncol(kmeans)-1)) {
  kmeans[kmeans[,i]>2.5,i]<-2.5
  kmeans[kmeans[,i]<-2.5,i]<- -2.5
}
##################################################
colnames(kmeans)<-gsub(colnames(kmeans),pattern = "_cibersort",replacement = "")
colnames(kmeans)<-gsub(colnames(kmeans),pattern = "_",replacement = " ")
##################################################
pp<-pheatmap(t(kmeans[,1:ncol(kmeans)-1]), annotation=annotation, 
             cluster_col=F,cluster_rows = TRUE,fontsize=12, 
             fontsize_row=12.5,show_colnames = F,cellheight = 13,
             treeheight_row=30,treeheight_col=30,width=2,heiht=20,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100)) 
#####################################################
#survival-analysis-of-TMEcluster
#####################################################
summary(km$cluster)
aa<-data.frame(km$cluster)
aa$ID<-rownames(aa)
aa$cluster<-as.factor(aa$km.cluster)
summary(aa$cluster)
###########################################
pdata_sur<-pdata[,c("ID","time","status")]
aa<-aa[aa$ID%in%pdata_sur$ID,];dim(aa)
pdata_sur<-pdata_sur[pdata_sur$ID%in%aa$ID,]
bb<-merge(pdata_sur,aa,by.x = "ID", by.y = "ID", all.x = T,all.y= T, incomparables = NULL)
head(bb)
summary(as.factor(bb$cluster))
############################################
bb$cluster<-revalue(bb$cluster,c("1"="TMEclusterB",
                                 "2"="TMEclusterA"))
bb<-bb[,c("ID","time","status","cluster")]
###########################################
bb<-bb[!is.na(bb$status),]
bb<-bb[bb$time>0,]
bb<-bb[!is.na(bb$cluster),]
#######################################
print(bb$cluster)
bb$cluster<-as.character(bb$cluster)
bb<-bb[order(bb$cluster,decreasing = T),]
y<-Surv(bb$time,bb$status)
summary(coxph(y~cluster,data = bb))
###########################################
fit<-survfit(y~bb$cluster,data =bb)
pp<-ggsurvplot(fit,data = bb,censor = TRUE, pval = T, 
               pval.size = 5,ncensor.plot = FALSE,
               conf.int = F,xlim = c(0:25),xlab = "Time in months", 
               break.time.by = 5,legend.labs = c('TMEclusterA','TMEclusterB'),
               submain="Cluster of TME cell in mUC",
               risk.table = T,palette = c("#2E9FDF","#E7B800","#E31A1C"));pp
#####################################
