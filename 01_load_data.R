library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gprofiler2)
library(stringr)
library(FlexDotPlot)
library(RColorBrewer)
library(ggpubr)
library(clusterProfiler)
library(DoubletFinder)
install.packages('hdf5r')
library(hdf5r)
aHC1.data<-Read10X_h5('/Users/apple/Desktop/study/fzxlab/scRNA/JIA/hlw_control/HC1_filtered_feature_bc_matrix.h5')
aHC2.data<-Read10X_h5('/Users/apple/Desktop/study/fzxlab/scRNA/JIA/hlw_control/HC2_filtered_feature_bc_matrix.h5')
aHC3.data<-Read10X_h5('/Users/apple/Desktop/study/fzxlab/scRNA/JIA/hlw_control/HC3_filtered_feature_bc_matrix.h5')
aHC4.data<-Read10X_h5('/Users/apple/Desktop/study/fzxlab/scRNA/JIA/hlw_control/HC4_filtered_feature_bc_matrix.h5')
aHC5.data<-Read10X_h5('/Users/apple/Desktop/study/fzxlab/scRNA/JIA/hlw_control/HC5_filtered_feature_bc_matrix.h5')

DHYP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/negative/DHY-P")
DYCP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/negative/DYC-P")
PBQZ.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/negative/PBQZ")
SZERP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/positive/JIA-JIA/SZER-P")
WJYP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/positive/JIA-JIA/WJY-P")
YJYP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/positive/JIA-JIA/YJY-S6")
ZJP.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/positive/JIA-JIA/ZJ-P")
CX.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/normal/CX")
LYR.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/normal/LYR")
WMJ.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/normal/WMJ")

SLE1.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/SLE/YE110-1")
SLE2.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/SLE/YE110-2")
SLE3.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/SLE/YE110-3")
SLE4.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/SLE/YE110-4")
pSS1.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/pSS/pSS_1")
pSS2.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/pSS/pSS_2")
pSS3.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/pSS/pSS_3")
pSS4.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/pSS/pSS_4")
pSS5.data<-Read10X(data.dir = "/Users/apple/Desktop/study/fzxlab/scRNA/JIA/第二篇查数据/pSS/pSS_5")

negative_name<-c('DHYP','DYCP','PBQZ')
positive_name<-c('SZERP','WJYP','YJYP','ZJP')
chc_name<-c('CX','LYR','WMJ')
ahc_name<-c('aHC1','aHC2','aHC3','aHC4','aHC5')
sle_name<-c('SLE1','SLE2','SLE3','SLE4')
pss_name<-c('pSS1','pSS2','pSS3','pSS4','pSS5')

negative_data<-list()
positive_data<-list()
chc_data<-list()
ahc_data<-list()
sle_data<-list()
pss_data<-list()

for (i in 1:length(negative_name)) {
  a<-paste0(negative_name[i],'.data')
  negative_data[[negative_name[i]]]<-expr_rename_10x(get(a),negative_name[i])
  negative_data[[i]]$meta<-umi_gene_count(negative_data[[i]]$expr)
  negative_data[[i]]$meta$group<-c(rep('b27 negative',length(negative_data[[i]]$meta$umi_count)))
  negative_data[[i]]$meta$person<-c(rep(negative_name[i],length(negative_data[[i]]$meta$umi_count)))
}

for (i in 1:length(positive_name)) {
  a<-paste0(positive_name[i],'.data')
  positive_data[[positive_name[i]]]<-expr_rename_10x(get(a),positive_name[i])
  positive_data[[i]]$meta<-umi_gene_count(positive_data[[i]]$expr)
  positive_data[[i]]$meta$group<-c(rep('b27 positive',length(positive_data[[i]]$meta$umi_count)))
  positive_data[[i]]$meta$person<-c(rep(positive_name[i],length(positive_data[[i]]$meta$umi_count)))
}

for (i in 1:length(chc_name)) {
  a<-paste0(chc_name[i],'.data')
  chc_data[[chc_name[i]]]<-expr_rename_10x(get(a),chc_name[i])
  chc_data[[i]]$meta<-umi_gene_count(chc_data[[i]]$expr)
  chc_data[[i]]$meta$group<-c(rep('cHC',length(chc_data[[i]]$meta$umi_count)))
  chc_data[[i]]$meta$person<-c(rep(chc_name[i],length(chc_data[[i]]$meta$umi_count)))
}

for (i in 1:length(ahc_name)) {
  a<-paste0(ahc_name[i],'.data')
  ahc_data[[ahc_name[i]]]<-expr_rename_10x(get(a),ahc_name[i])
  ahc_data[[i]]$meta<-umi_gene_count(ahc_data[[i]]$expr)
  ahc_data[[i]]$meta$group<-c(rep('aHC',length(ahc_data[[i]]$meta$umi_count)))
  ahc_data[[i]]$meta$person<-c(rep(ahc_name[i],length(ahc_data[[i]]$meta$umi_count)))
}

for (i in 1:length(sle_name)) {
  a<-paste0(sle_name[i],'.data')
  sle_data[[sle_name[i]]]<-expr_rename_10x(get(a),sle_name[i])
  sle_data[[i]]$meta<-umi_gene_count(sle_data[[i]]$expr)
  sle_data[[i]]$meta$group<-c(rep('SLE',length(sle_data[[i]]$meta$umi_count)))
  sle_data[[i]]$meta$person<-c(rep(sle_name[i],length(sle_data[[i]]$meta$umi_count)))
}

for (i in 1:length(pss_name)) {
  a<-paste0(pss_name[i],'.data')
  pss_data[[pss_name[i]]]<-expr_rename_10x(get(a),pss_name[i])
  pss_data[[i]]$meta<-umi_gene_count(pss_data[[i]]$expr)
  pss_data[[i]]$meta$group<-c(rep('pSS',length(pss_data[[i]]$meta$umi_count)))
  pss_data[[i]]$meta$person<-c(rep(pss_name[i],length(pss_data[[i]]$meta$umi_count)))
}

all_expr<-cbind(negative_data$DHYP$expr,negative_data$DYCP$expr,negative_data$PBQZ$expr,
                positive_data$SZERP$expr,positive_data$WJYP$expr,positive_data$YJYP$expr,positive_data$ZJP$expr,
                chc_data$CX$expr,chc_data$LYR$expr,chc_data$WMJ$expr,
                ahc_data$aHC1$expr,ahc_data$aHC2$expr,ahc_data$aHC3$expr,ahc_data$aHC4$expr,ahc_data$aHC5$expr,
                pss_data$pSS1$expr,pss_data$pSS2$expr,pss_data$pSS3$expr,pss_data$pSS4$expr,pss_data$pSS5$expr,
                sle_data$SLE1$expr,sle_data$SLE2$expr,sle_data$SLE3$expr,sle_data$SLE4$expr)
all_meta<-rbind(negative_data$DHYP$meta,negative_data$DYCP$meta,negative_data$PBQZ$meta,
                positive_data$SZERP$meta,positive_data$WJYP$meta,positive_data$YJYP$meta,positive_data$ZJP$meta,
                chc_data$CX$meta,chc_data$LYR$meta,chc_data$WMJ$meta,
                ahc_data$aHC1$meta,ahc_data$aHC2$meta,ahc_data$aHC3$meta,ahc_data$aHC4$meta,ahc_data$aHC5$meta,
                pss_data$pSS1$meta,pss_data$pSS2$meta,pss_data$pSS3$meta,pss_data$pSS4$meta,pss_data$pSS5$meta,
                sle_data$SLE1$meta,sle_data$SLE2$meta,sle_data$SLE3$meta,sle_data$SLE4$meta)
