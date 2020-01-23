




#ource-packages
library(IOBR)
library(tidyverse)
library(GSVA)
##########################

#load expression profile
(load(paste0("log2TPM-symbol-IMvigor210.RData")))
###########################

#calculate prevelent signatures score(see Supplementary Table S11)
#calculate_sig_score_zdq a function in an unplished R package (IOBR) developed by Dongqiang Zeng
pdata_sig<-calculate_sig_score_zdq(pdata = pdata,
                                 eset  = eset,
                                 signature = my_signatures,
                                 method = "pca",
                                 mingenecounts = 2)
##########################

#calculate signatures score of gene sets obtained from MsigDB: (http://software.broadinstitute.org/gsea/msigdb/index.jsp)
pdata_sig<-calculate_sig_score_zdq(pdata = pdata,
                                   eset = eset,
                                   signature =c(hallmark_gene_set,
                                                go_bp_gene_set,
                                                go_cc_gene_set,
                                                go_mf_gene_set,
                                                kegg_gene_set,
                                                reactome_gene_set,
                                                cm_gene_set),
                                   method = "ssgsea",
                                   mingenecounts = 2)

###########################






