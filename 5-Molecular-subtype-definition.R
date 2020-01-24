

#loading R packages
##################################
# devtools::install_github("cit-bioinfo/consensusMIBC", build_vignettes = TRUE)
##################################
library(consensusMIBC)
# Single sample classification
##########################
(load(paste0("5-TPM_entrezid-IMvigor210.RData")))
eset<-log2(expMatrix_tpm+1)
eset<-as.matrix(eset[order(rownames(eset),decreasing = F),])
########################
res <- getConsensusClass(eset,minCor = 0.2,gene_id = "entrezgene")
head(res)
table(res$consensusClass)
res<-rownames_to_column(res,var = "ID")
##########################
# devtools::install_github("cit-bioinfo/BLCAsubtyping", build_vignettes = TRUE)
library(BLCAsubtyping)
vignette("BLCAsubtyping")
data(cit) 
#########################
(load(paste0(path_uc,"6-TPM-symbol-IMvigor210.RData")))
eset<-log2(eset+1)
eset<-as.matrix(eset[order(rownames(eset),decreasing = F),])
eset[1:7,1:7]
summary(duplicated(rownames(eset)))
sum(is.na(eset))
max(eset);min(eset)
summary(rownames(eset)%in%rownames(cit))
summary(rownames(cit)%in%rownames(eset))
eset<-eset[rownames(eset)%in%rownames(cit),]
########################
cl <- classify(expMat = as.matrix(eset), 
               symbol = "Gene.Symbol",
               classification.systems = c("Baylor", "UNC", "MDA", "CIT", "Lund", "TCGA"))
head(cl)
########################
subtype<-merge(cl,res,by="ID",all = F)
head(subtype)
#######################