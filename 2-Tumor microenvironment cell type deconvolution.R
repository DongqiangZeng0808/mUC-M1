

#source packages
library(IOBR)
library(xCell)
library(EPIC)
library(MCPcounter)
library(tidyverse)
library(limSolve)
##########################

#load expression profile
(load(paste0("log2TPM-symbol-IMvigor210.RData")))
###########################
eset<-as.data.frame(eset)
eset<-rownames_to_column(eset,var = "symbol")
###########################

#CIBERSORT
result<-IOBR::deconvolute(gene_expression = eset,method = "cibersort",tumor = T,arrays = array_judge)
result<-column_to_rownames(result,var = "cell_type")
result<-as.data.frame(t(result))
colnames(result)<-paste(colnames(result),"_cibersort",sep = "")
cibersort<-rownames_to_column(result,var = "ID")
############################
############################
eset1<-column_to_rownames(eset,var = "symbol")
eset1<-eset1[!duplicated(rownames(eset1)),]
eset1<-eset1[!is.na(rownames(eset1)),]
#############################

#TIMER
##############################
result<-IOBR::deconvolute(gene_expression = eset1,method = "timer",indications = rep(tumor_type,dim(eset1)[2]))
head(result);dim(result)
rownames(result)<-NULL
result<-column_to_rownames(result,var = "cell_type")
result<-as.data.frame(t(result))
colnames(result)<-paste(colnames(result),"_TIMER",sep = "")
timer<-rownames_to_column(result,var = "ID")
#############################

#xCell
##############################
xcellresult<-IOBR::xcell_deconvo_zdq(ProjectID = ProjectID, eset = eset1 )
patien_number<-dim(xcellresult)[1]
xcellresult[1:10,1:10]
####################################

#MCP-counter
####################################
mcpresult<-IOBR::mcp_deconvo_zdq(ProjectID = ProjectID, eset = eset1 )
patien_number<-dim(mcpresult)[1]
mcpresult[1:10,1:10]
####################################

#EPIC
####################################
epicresult<-IOBR::epic_deconvo_zdq(ProjectID = ProjectID, eset = eset1 )
patien_number<-dim(epicresult)[1]
epicresult[1:10,1:10]
####################################










