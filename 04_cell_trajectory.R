setwd('D:\\scRNA\\JIA\\trajectory')
library(monocle)
library(pheatmap)
library(scales)
library(ggplotify)
library(viridis)
library(ggpubr)
library(Seurat)
devtools::load_all('E:\\R-4.2.1\\library\\monocle')

saveRDS(subset(seurat,subset = group_2 == 'aHC' & celltype_1 == 'T cells'),
        'D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_t.rds')
saveRDS(subset(seurat,subset = group_2 == 'aHC' & celltype_1 == 'B cells'),
        'D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_b.rds')
saveRDS(subset(seurat,subset = group_2 == 'aHC' & celltype_1 == 'Myeloid cells'),
        'D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_m.rds')
saveRDS(subset(seurat,subset = group_2 == 'cHC' & celltype_1 == 'T cells'),
        'D:\\scRNA\\JIA\\trajectory\\cHC\\chc_t.rds')
saveRDS(subset(seurat,subset = group_2 == 'cHC' & celltype_1 == 'B cells'),
        'D:\\scRNA\\JIA\\trajectory\\cHC\\chc_b.rds')
saveRDS(subset(seurat,subset = group_2 == 'cHC' & celltype_1 == 'Myeloid cells'),
        'D:\\scRNA\\JIA\\trajectory\\cHC\\chc_m.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27-' & celltype_1 == 'T cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_t.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27-' & celltype_1 == 'B cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_b.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27-' & celltype_1 == 'Myeloid cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_m.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27+' & celltype_1 == 'T cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_t.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27+' & celltype_1 == 'B cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_b.rds')
saveRDS(subset(seurat,subset = group == 'HLA-B27+' & celltype_1 == 'Myeloid cells'),
        'D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_m.rds')

main_figure[['f2_a1']]<-DimPlot(subset(jia_seurat,subset = celltype_1 == 'T cells'),group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(6,'Set2')))+
  labs(title = 'JIA & cHC PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure$f2_a1

jia_t<-subset(jia_t,subset = group_2 == 'JIA')
jia_t<-NormalizeData(jia_t)
jia_t<-FindVariableFeatures(jia_t,selection.method = 'vst',nfeatures = 2000)
jia_t<-ScaleData(jia_t)
jia_t<-RunPCA(jia_t,features = VariableFeatures(object = jia_t))
ElbowPlot(jia_t)
jia_t<-FindNeighbors(jia_t,dims = 1:10)
jia_t<-FindClusters(jia_t,resolution = 0.3)
jia_t<-RunUMAP(jia_t,dims = 1:10)
main_figure[['f2_a1']]<-DimPlot(jia_t,group.by = 'celltype_2',label = TRUE,cols = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  labs(title = 'JIA PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_a1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_A1.pdf',units = 'in',dpi = 300,
       width = 7,height = 7,device = "pdf")

jia_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\JIA\\jia_t.rds')
jia_t_expr_matrix<-as(as.matrix(subset(jia_t,subset = group_2 == 'JIA')@assays$RNA@counts),'sparseMatrix')
jia_t_p_data<-subset(jia_t,subset = group_2 == 'JIA')@meta.data
jia_t_f_data<-data.frame(gene_short_name=row.names(subset(jia_t,subset = group_2 == 'JIA')),
                         row.names = row.names(subset(jia_t,subset = group_2 == 'JIA')))
jia_t_pd<-new('AnnotatedDataFrame',data=jia_t_p_data)
jia_t_fd<-new('AnnotatedDataFrame',data=jia_t_f_data)
jia_t_cds<-newCellDataSet(jia_t_expr_matrix,
                          phenoData = jia_t_pd,
                          featureData = jia_t_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
jia_t_cds<-estimateSizeFactors(jia_t_cds)
jia_t_cds<-estimateDispersions(jia_t_cds)
jia_t_cds<-detectGenes(jia_t_cds,min_expr = 0.5)
disp_table<-dispersionTable(jia_t_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_t_cds<-setOrderingFilter(jia_t_cds,ordering_genes)
jia_t_cds<-reduceDimension(jia_t_cds,max_components = 2,method = 'DDRTree')
jia_t_cds<-orderCells(jia_t_cds,reverse = F)
main_figure<-list()
main_figure[['f2_b']]<-plot_cell_trajectory(jia_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  #scale_fill_manual(values = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  ggtitle('JIA PBMC T cell subtypes trajectory')
main_figure[['f2_b']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F2_C1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_b2']]<-plot_cell_trajectory(jia_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('JIA PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f2_b2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F2_B2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

jia_t_disp_table<-dispersionTable(jia_t_cds)
jia_t_ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_t_time_diff<-differentialGeneTest(jia_t_cds[jia_t_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
jia_t_time_diff<-jia_t_time_diff[,c(5,2,3,4,1,6,7)]
#jia_t_time_genes<-top_n(jia_t_time_diff, n=5, desc(qval)) %>% pull(gene_short_name) %>% as.character()
#jia_t_time_genes<-arrange(jia_t_time_diff,jia_t_time_diff$qval)
#jia_t_time_genes<-row.names(jia_t_time_genes[1:20,])

jia_t_markers<-FindAllMarkers(subset(jia_t,subset = group_2 == 'JIA'),only.pos = TRUE,logfc.threshold = 0.5)
jia_t_top10<-jia_t_markers %>% group_by(cluster) %>% top_n(n=7,wt=avg_log2FC)
jia_t_top10_ordergene<-jia_t_time_diff[jia_t_top10$gene,]
jia_t_time_genes<-jia_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
jia_t_time_genes<-unique(jia_t_time_genes)
jia_t_time_genes<-na.omit(jia_t_time_genes)
# plot_pseudotime_heatmap(jia_t_cds[jia_t_time_genes,],num_clusters = 1,show_rownames = T,
#                         return_heatmap = T,cluster_rows = FALSE,
#                         hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

jia_t_newdata<-data.frame(Pseudotime = seq(min(pData(jia_t_cds)$Pseudotime),
                                           max(pData(jia_t_cds)$Pseudotime),length.out=100))
jia_t_m<-genSmoothCurves(jia_t_cds[jia_t_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = jia_t_newdata)
jia_t_m<-jia_t_m[!apply(jia_t_m,1,sum)==0,]
jia_t_m<-log10(jia_t_m+1)
jia_t_m<-jia_t_m[!apply(jia_t_m,1,sd)==0,]
jia_t_m<-Matrix::t(scale(Matrix::t(jia_t_m),center = TRUE))
jia_t_m<-jia_t_m[is.na(row.names(jia_t_m))==FALSE,]
jia_t_m[is.nan(jia_t_m)]=0
jia_t_m[jia_t_m>3]=3
jia_t_m[jia_t_m<-3]=-3
jia_t_row_dist<-as.dist((1-cor(Matrix::t(jia_t_m)))/2)
jia_t_row_dist[is.na(jia_t_row_dist)]<-1
p1<-pheatmap(jia_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = jia_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
jia_t_annotation_col<-data.frame(pseudotime=rescale(jia_t_newdata$Pseudotime,to=c(-1,1)))
row.names(jia_t_annotation_col)<-colnames(jia_t_m)
jia_t_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(jia_t_annotation_row)<-rownames(jia_t_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

jia_t_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f2_c']]<-pheatmap(jia_t_m,useRaster = T,cluster_rows = TRUE,
                                cluster_cols = FALSE,show_rownames = TRUE,
                                show_colnames = FALSE,cutree_rows = 4,
                                clustering_distance_rows = jia_t_row_dist,
                                border_color = NA,filename = NA,
                                color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                annotation_col = jia_t_annotation_col,
                                annotation_colors = jia_t_anno_colors,
                                annotation_row = jia_t_annotation_row)
main_figure[['f2_c']]<-as.ggplot(main_figure[['f2_c']])
main_figure[['f2_c']]<-ggarrange(main_figure[['f2_c']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'left')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f2_c']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_C.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

chc_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\cHC\\chc_t.rds')
chc_t<-NormalizeData(chc_t)
chc_t<-FindVariableFeatures(chc_t,selection.method = 'vst',nfeatures = 2000)
chc_t<-ScaleData(chc_t)
chc_t<-RunPCA(chc_t,features = VariableFeatures(object = chc_t))
ElbowPlot(chc_t)
chc_t<-FindNeighbors(chc_t,dims = 1:10)
chc_t<-FindClusters(chc_t,resolution = 0.3)
chc_t<-RunUMAP(chc_t,dims = 1:10)
main_figure[['f2_b1']]<-DimPlot(chc_t,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(6,'Set1')))+
  labs(title = 'cHC PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_b1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F2_B1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

chc_t_expr_matrix<-as(as.matrix(chc_t@assays$RNA@counts),'sparseMatrix')
chc_t_p_data<-chc_t@meta.data
chc_t_f_data<-data.frame(gene_short_name=row.names(chc_t),row.names = row.names(chc_t))
chc_t_pd<-new('AnnotatedDataFrame',data=chc_t_p_data)
chc_t_fd<-new('AnnotatedDataFrame',data=chc_t_f_data)
chc_t_cds<-newCellDataSet(chc_t_expr_matrix,
                          phenoData = chc_t_pd,
                          featureData = chc_t_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
chc_t_cds<-estimateSizeFactors(chc_t_cds)
chc_t_cds<-estimateDispersions(chc_t_cds)
chc_t_cds<-detectGenes(chc_t_cds,min_expr = 0.5)
disp_table<-dispersionTable(chc_t_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_t_cds<-setOrderingFilter(chc_t_cds,ordering_genes)
chc_t_cds<-reduceDimension(chc_t_cds,max_components = 2,method = 'DDRTree',
                           num_dim = 6)
chc_t_cds<-orderCells(chc_t_cds,reverse = F)
#main_figure<-list()
main_figure[['f2_b2']]<-plot_cell_trajectory(chc_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  ggtitle('cHC PBMC T cell subtypes trajectory')
main_figure[['f2_b2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F2_E1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_b3']]<-plot_cell_trajectory(chc_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('cHC PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f2_b3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_B3.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

chc_t_disp_table<-dispersionTable(chc_t_cds)
chc_t_ordering_genes<-as.character(subset(chc_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_t_time_diff<-differentialGeneTest(chc_t_cds[chc_t_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
chc_t_time_diff<-chc_t_time_diff[,c(5,2,3,4,1,6)]
chc_t_markers<-FindAllMarkers(chc_t,only.pos = TRUE,logfc.threshold = 0.5)
chc_t_top10<-chc_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
chc_t_top10_ordergene<-chc_t_time_diff[chc_t_top10$gene,]
chc_t_time_genes<-chc_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
chc_t_time_genes<-unique(chc_t_time_genes)
chc_t_time_genes<-na.omit(chc_t_time_genes)
chc_t_newdata<-data.frame(Pseudotime = seq(min(pData(chc_t_cds)$Pseudotime),
                                           max(pData(chc_t_cds)$Pseudotime),length.out=100))
chc_t_m<-genSmoothCurves(chc_t_cds[chc_t_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = chc_t_newdata)
chc_t_m<-chc_t_m[!apply(chc_t_m,1,sum)==0,]
chc_t_m<-log10(chc_t_m+1)
chc_t_m<-chc_t_m[!apply(chc_t_m,1,sd)==0,]
chc_t_m<-Matrix::t(scale(Matrix::t(chc_t_m),center = TRUE))
chc_t_m<-chc_t_m[is.na(row.names(chc_t_m))==FALSE,]
chc_t_m[is.nan(chc_t_m)]=0
chc_t_m[chc_t_m>3]=3
chc_t_m[chc_t_m<-3]=-3
chc_t_row_dist<-as.dist((1-cor(Matrix::t(chc_t_m)))/2)
chc_t_row_dist[is.na(chc_t_row_dist)]<-1
p1<-pheatmap(chc_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = chc_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
chc_t_annotation_col<-data.frame(pseudotime=rescale(chc_t_newdata$Pseudotime,to=c(-1,1)))
row.names(chc_t_annotation_col)<-colnames(chc_t_m)
chc_t_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(chc_t_annotation_row)<-rownames(chc_t_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

chc_t_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['chc_t']]<-pheatmap(chc_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = chc_t_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = chc_t_annotation_col,
                                 annotation_colors = chc_t_anno_colors,
                                 annotation_row = chc_t_annotation_row)
main_figure[['chc_t']]<-as.ggplot(main_figure[['chc_t']])
main_figure[['chc_t']]<-ggarrange(main_figure[['chc_t']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['chc_t']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\chc_t_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

jia_b<-subset(jia_b,subset = group_2 == 'JIA')
jia_b<-NormalizeData(jia_b)
jia_b<-FindVariableFeatures(jia_b,selection.method = 'vst',nfeatures = 2000)
jia_b<-ScaleData(jia_b)
jia_b<-RunPCA(jia_b,features = VariableFeatures(object = jia_b))
ElbowPlot(jia_b)
jia_b<-FindNeighbors(jia_b,dims = 1:10)
jia_b<-FindClusters(jia_b,resolution = 0.3)
jia_b<-RunUMAP(jia_b,dims = 1:10)
main_figure[['f2_d1']]<-DimPlot(jia_b,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(3,'Set1')))+
  labs(title = 'JIA PBMC B cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_d1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_D1.pdf',units = 'in',dpi = 300,
       width = 7,height = 7,device = "pdf")

jia_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\JIA\\jia_b.rds')
jia_b_expr_matrix<-as(as.matrix(subset(jia_b,subset = group_2 == 'JIA')@assays$RNA@counts),'sparseMatrix')
jia_b_p_data<-subset(jia_b,subset = group_2 == 'JIA')@meta.data
jia_b_f_data<-data.frame(gene_short_name=row.names(subset(jia_b,subset = group_2 == 'JIA')),
                         row.names = row.names(subset(jia_b,subset = group_2 == 'JIA')))
jia_b_pd<-new('AnnotatedDataFrame',data=jia_b_p_data)
jia_b_fd<-new('AnnotatedDataFrame',data=jia_b_f_data)
jia_b_cds<-newCellDataSet(jia_b_expr_matrix,
                          phenoData = jia_b_pd,
                          featureData = jia_b_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
jia_b_cds<-estimateSizeFactors(jia_b_cds)
jia_b_cds<-estimateDispersions(jia_b_cds)
jia_b_cds<-detectGenes(jia_b_cds,min_expr = 0.5)
disp_table<-dispersionTable(jia_b_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_b_cds<-setOrderingFilter(jia_b_cds,ordering_genes)
jia_b_cds<-reduceDimension(jia_b_cds,max_components = 2,method = 'DDRTree')
jia_b_cds<-orderCells(jia_b_cds,reverse = F)
main_figure[['f2_e']]<-plot_cell_trajectory(jia_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#EC8D63','#DE4B3F','#DE4247'))+
  ggtitle('JIA PBMC B cell subtypes trajectory')
main_figure[['f2_e']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_C1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_e2']]<-plot_cell_trajectory(jia_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('JIA PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f2_e2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-10-12\\F5_C2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

jia_b_disp_table<-dispersionTable(jia_b_cds)
jia_b_ordering_genes<-as.character(subset(jia_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_b_time_diff<-differentialGeneTest(jia_b_cds[jia_b_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
jia_b_time_diff<-jia_b_time_diff[,c(5,2,3,4,1,6,7)]
jia_b_markers<-FindAllMarkers(subset(jia_b,subset = group_2 == 'JIA'),only.pos = TRUE,logfc.threshold = 0.5)
jia_b_top10<-jia_b_markers %>% group_by(cluster) %>% top_n(n=8,wt=avg_log2FC)
jia_b_top10_ordergene<-jia_b_time_diff[jia_b_top10$gene,]
jia_b_time_genes<-jia_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
jia_b_time_genes<-unique(jia_b_time_genes)
jia_b_time_genes<-na.omit(jia_b_time_genes)
jia_b_newdata<-data.frame(Pseudotime = seq(min(pData(jia_b_cds)$Pseudotime),
                                           max(pData(jia_b_cds)$Pseudotime),length.out=100))
jia_b_m<-genSmoothCurves(jia_b_cds[jia_b_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = jia_b_newdata)
jia_b_m<-jia_b_m[!apply(jia_b_m,1,sum)==0,]
jia_b_m<-log10(jia_b_m+1)
jia_b_m<-jia_b_m[!apply(jia_b_m,1,sd)==0,]
jia_b_m<-Matrix::t(scale(Matrix::t(jia_b_m),center = TRUE))
jia_b_m<-jia_b_m[is.na(row.names(jia_b_m))==FALSE,]
jia_b_m[is.nan(jia_b_m)]=0
jia_b_m[jia_b_m>3]=3
jia_b_m[jia_b_m<-3]=-3
jia_b_row_dist<-as.dist((1-cor(Matrix::t(jia_b_m)))/2)
jia_b_row_dist[is.na(jia_b_row_dist)]<-1

p1<-pheatmap(jia_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = jia_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1

jia_b_annotation_col<-data.frame(pseudotime=rescale(jia_b_newdata$Pseudotime,to=c(-1,1)))
row.names(jia_b_annotation_col)<-colnames(jia_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(jia_b_m)
jia_b_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f2_f']]<-pheatmap(jia_b_m,useRaster = T,cluster_rows = TRUE,
                                cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                show_colnames = FALSE,cutree_rows = 4,
                                clustering_distance_rows = jia_b_row_dist,
                                border_color = NA,filename = NA,
                                color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                annotation_col = jia_b_annotation_col,
                                annotation_colors = jia_b_anno_colors,
                                annotation_row = annotation_row)
main_figure[['f2_f']]<-as.ggplot(main_figure[['f2_f']])
main_figure[['f2_f']]<-ggarrange(main_figure[['f2_f']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f2_f']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_F.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

chc_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\cHC\\chc_b.rds')
chc_b<-NormalizeData(chc_b)
chc_b<-FindVariableFeatures(chc_b,selection.method = 'vst',nfeatures = 2000)
chc_b<-ScaleData(chc_b)
chc_b<-RunPCA(chc_b,features = VariableFeatures(object = chc_b))
ElbowPlot(chc_b)
chc_b<-FindNeighbors(chc_b,dims = 1:20)
chc_b<-FindClusters(chc_b,resolution = 0.3)
chc_b<-RunUMAP(chc_b,dims = 1:20)
main_figure[['f2_e1']]<-DimPlot(chc_b,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(3,'Set1')))+
  labs(title = 'cHC PBMC B cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_e1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_E1.pdf',units = 'in',dpi = 300,
       width = 7,height = 7,device = "pdf")

chc_b_expr_matrix<-as(as.matrix(chc_b@assays$RNA@counts),'sparseMatrix')
chc_b_p_data<-chc_b@meta.data
chc_b_f_data<-data.frame(gene_short_name=row.names(chc_b),row.names = row.names(chc_b))
chc_b_pd<-new('AnnotatedDataFrame',data=chc_b_p_data)
chc_b_fd<-new('AnnotatedDataFrame',data=chc_b_f_data)
chc_b_cds<-newCellDataSet(chc_b_expr_matrix,
                          phenoData = chc_b_pd,
                          featureData = chc_b_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
chc_b_cds<-estimateSizeFactors(chc_b_cds)
chc_b_cds<-estimateDispersions(chc_b_cds)
chc_b_cds<-detectGenes(chc_b_cds,min_expr = 0.5)
disp_table<-dispersionTable(chc_b_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_b_cds<-setOrderingFilter(chc_b_cds,ordering_genes)
chc_b_cds<-reduceDimension(chc_b_cds,max_components = 2,method = 'DDRTree',num_dim = 6)
chc_b_cds<-orderCells(chc_b_cds,reverse = F)
main_figure[['f2_e2']]<-plot_cell_trajectory(chc_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#EC8D63','#DE4B3F','#DE4247'))+
  ggtitle('cHC PBMC B cell subtypes trajectory')
main_figure[['f2_e2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_E1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_e3']]<-plot_cell_trajectory(chc_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('cHC PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f2_e3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_E2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

chc_b_disp_table<-dispersionTable(chc_b_cds)
chc_b_ordering_genes<-as.character(subset(chc_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_b_time_diff<-differentialGeneTest(chc_b_cds[chc_b_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
chc_b_time_diff<-chc_b_time_diff[,c(5,2,3,4,1,6,7)]
chc_b_markers<-FindAllMarkers(chc_b,only.pos = TRUE,logfc.threshold = 0.5)
chc_b_top10<-chc_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
chc_b_top10_ordergene<-chc_b_time_diff[chc_b_top10$gene,]
chc_b_time_genes<-chc_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
chc_b_time_genes<-unique(chc_b_time_genes)
chc_b_time_genes<-na.omit(chc_b_time_genes)
chc_b_newdata<-data.frame(Pseudotime = seq(min(pData(chc_b_cds)$Pseudotime),
                                           max(pData(chc_b_cds)$Pseudotime),length.out=100))
chc_b_m<-genSmoothCurves(chc_b_cds[chc_b_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = chc_b_newdata)
chc_b_m<-chc_b_m[!apply(chc_b_m,1,sum)==0,]
chc_b_m<-log10(chc_b_m+1)
chc_b_m<-chc_b_m[!apply(chc_b_m,1,sd)==0,]
chc_b_m<-Matrix::t(scale(Matrix::t(chc_b_m),center = TRUE))
chc_b_m<-chc_b_m[is.na(row.names(chc_b_m))==FALSE,]
chc_b_m[is.nan(chc_b_m)]=0
chc_b_m[chc_b_m>3]=3
chc_b_m[chc_b_m<-3]=-3
chc_b_row_dist<-as.dist((1-cor(Matrix::t(chc_b_m)))/2)
chc_b_row_dist[is.na(chc_b_row_dist)]<-1
p1<-pheatmap(chc_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = chc_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
chc_b_annotation_col<-data.frame(pseudotime=rescale(chc_b_newdata$Pseudotime,to=c(-1,1)))
row.names(chc_b_annotation_col)<-colnames(chc_b_m)
chc_b_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(chc_b_annotation_row)<-rownames(chc_b_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

chc_b_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['chc_b']]<-pheatmap(chc_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = chc_b_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = chc_b_annotation_col,
                                 annotation_colors = chc_b_anno_colors,
                                 annotation_row = chc_b_annotation_row)
main_figure[['chc_b']]<-as.ggplot(main_figure[['chc_b']])
main_figure[['chc_b']]<-ggarrange(main_figure[['chc_b']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['chc_b']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\chc_b_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

jia_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\JIA\\jia_m.rds')
#jia_m<-subset(jia_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
jia_m<-subset(jia_m,subset = group_2 == 'JIA')
#jia_m<-subset(jia_seurat,subset = celltype_1 == 'Myeloid cells')
jia_m<-NormalizeData(jia_m)
jia_m<-FindVariableFeatures(jia_m,selection.method = 'vst',nfeatures = 2000)
jia_m<-ScaleData(jia_m)
jia_m<-RunPCA(jia_m,features = VariableFeatures(object = jia_m))
ElbowPlot(jia_m)
jia_m<-FindNeighbors(jia_m,dims = 1:10)
jia_m<-FindClusters(jia_m,resolution = 0.3)
jia_m<-RunUMAP(jia_m,dims = 1:10)
main_figure[['f2_g1']]<-DimPlot(jia_m,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(4,'Set1')))+
  labs(title = 'JIA PBMC Myeloid cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_g1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_G1.pdf',units = 'in',dpi = 300,
       width = 7,height = 7,device = "pdf")

jia_m_expr_matrix<-as(as.matrix(subset(jia_m,subset = group_2 == 'JIA')@assays$RNA@counts),'sparseMatrix')
jia_m_p_data<-subset(jia_m,subset = group_2 == 'JIA')@meta.data
jia_m_f_data<-data.frame(gene_short_name=row.names(subset(jia_m,subset = group_2 == 'JIA')),
                         row.names = row.names(subset(jia_m,subset = group_2 == 'JIA')))
jia_m_pd<-new('AnnotatedDataFrame',data=jia_m_p_data)
jia_m_fd<-new('AnnotatedDataFrame',data=jia_m_f_data)
jia_m_cds<-newCellDataSet(jia_m_expr_matrix,
                          phenoData = jia_m_pd,
                          featureData = jia_m_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
jia_m_cds<-estimateSizeFactors(jia_m_cds)
jia_m_cds<-estimateDispersions(jia_m_cds)
jia_m_cds<-detectGenes(jia_m_cds,min_expr = 0.5)
disp_table<-dispersionTable(jia_m_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_m_cds<-setOrderingFilter(jia_m_cds,ordering_genes)
jia_m_cds<-reduceDimension(jia_m_cds,max_components = 2,method = 'DDRTree')
jia_m_cds<-orderCells(jia_m_cds,reverse = F)
main_figure[['f2_h']]<-plot_cell_trajectory(jia_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#65A644','#657A51','#8DC591'))+
  ggtitle('JIA PBMC Monocyte subtypes trajectory')
main_figure[['f2_h']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F6_C1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_h2']]<-plot_cell_trajectory(jia_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('JIA PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f2_h2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F6_C2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

jia_m_disp_table<-dispersionTable(jia_m_cds)
jia_m_ordering_genes<-as.character(subset(jia_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
jia_m_time_diff<-differentialGeneTest(jia_m_cds[jia_m_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
jia_m_time_diff<-jia_m_time_diff[,c(5,2,3,4,1,6,7)]
jia_m_markers<-FindAllMarkers(subset(jia_m,subset = group_2 == 'JIA'),only.pos = TRUE,logfc.threshold = 0.5)
jia_m_top10<-jia_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
jia_m_top10_ordergene<-jia_m_time_diff[jia_m_top10$gene,]
jia_m_time_genes<-jia_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
jia_m_time_genes<-unique(jia_m_time_genes)
jia_m_time_genes<-na.omit(jia_m_time_genes)
jia_m_newdata<-data.frame(Pseudotime = seq(min(pData(jia_m_cds)$Pseudotime),
                                           max(pData(jia_m_cds)$Pseudotime),length.out=100))
jia_m_m<-genSmoothCurves(jia_m_cds[jia_m_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = jia_m_newdata)
jia_m_m<-jia_m_m[!apply(jia_m_m,1,sum)==0,]
jia_m_m<-log10(jia_m_m+1)
jia_m_m<-jia_m_m[!apply(jia_m_m,1,sd)==0,]
jia_m_m<-Matrix::t(scale(Matrix::t(jia_m_m),center = TRUE))
jia_m_m<-jia_m_m[is.na(row.names(jia_m_m))==FALSE,]
jia_m_m[is.nan(jia_m_m)]=0
jia_m_m[jia_m_m>3]=3
jia_m_m[jia_m_m<-3]=-3
jia_m_row_dist<-as.dist((1-cor(Matrix::t(jia_m_m)))/2)
jia_m_row_dist[is.na(jia_m_row_dist)]<-1

p1<-pheatmap(jia_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = jia_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

jia_m_annotation_col<-data.frame(pseudotime=rescale(jia_m_newdata$Pseudotime,to=c(-1,1)))
row.names(jia_m_annotation_col)<-colnames(jia_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(jia_m_m)
jia_m_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f2_i']]<-pheatmap(jia_m_m,useRaster = T,cluster_rows = TRUE,
                                cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                show_colnames = FALSE,cutree_rows = 4,
                                clustering_distance_rows = jia_m_row_dist,
                                border_color = NA,filename = NA,
                                color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                annotation_col = jia_m_annotation_col,
                                annotation_colors = jia_m_anno_colors,
                                annotation_row = annotation_row)
main_figure[['f2_i']]<-as.ggplot(main_figure[['f2_i']])
main_figure[['f2_i']]<-ggarrange(main_figure[['f2_i']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f2_i']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_I.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

chc_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\cHC\\chc_m.rds')
chc_m<-NormalizeData(chc_m)
chc_m<-FindVariableFeatures(chc_m,selection.method = 'vst',nfeatures = 2000)
chc_m<-ScaleData(chc_m)
chc_m<-RunPCA(chc_m,features = VariableFeatures(object = chc_m))
ElbowPlot(chc_m)
chc_m<-FindNeighbors(chc_m,dims = 1:10)
chc_m<-FindClusters(chc_m,resolution = 0.3)
chc_m<-RunUMAP(chc_m,dims = 1:10)
main_figure[['f2_h1']]<-DimPlot(chc_m,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(4,'Set1')))+
  labs(title = 'cHC PBMC Myeloid cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f2_h1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_H1.pdf',units = 'in',dpi = 300,
       width = 7,height = 7,device = "pdf")

chc_m<-subset(chc_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
chc_m_expr_matrix<-as(as.matrix(chc_m@assays$RNA@counts),'sparseMatrix')
chc_m_p_data<-chc_m@meta.data
chc_m_f_data<-data.frame(gene_short_name=row.names(chc_m),row.names = row.names(chc_m))
chc_m_pd<-new('AnnotatedDataFrame',data=chc_m_p_data)
chc_m_fd<-new('AnnotatedDataFrame',data=chc_m_f_data)
chc_m_cds<-newCellDataSet(chc_m_expr_matrix,
                          phenoData = chc_m_pd,
                          featureData = chc_m_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
chc_m_cds<-estimateSizeFactors(chc_m_cds)
chc_m_cds<-estimateDispersions(chc_m_cds)
chc_m_cds<-detectGenes(chc_m_cds,min_expr = 0.5)
disp_table<-dispersionTable(chc_m_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_m_cds<-setOrderingFilter(chc_m_cds,ordering_genes)
chc_m_cds<-reduceDimension(chc_m_cds,max_components = 2,method = 'DDRTree',num_dim = 6)
chc_m_cds<-orderCells(chc_m_cds,reverse = F)
main_figure[['f2_h2']]<-plot_cell_trajectory(chc_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#65A644','#657A51','#8DC591'))+
  ggtitle('cHC PBMC Monocyte subtypes trajectory')
main_figure[['f2_h2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F6_E1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f2_h3']]<-plot_cell_trajectory(chc_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('cHC PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f2_h3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F2_H3.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

chc_m_disp_table<-dispersionTable(chc_m_cds)
chc_m_ordering_genes<-as.character(subset(chc_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
chc_m_time_diff<-differentialGeneTest(chc_m_cds[chc_m_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
chc_m_time_diff<-chc_m_time_diff[,c(5,2,3,4,1,6,7)]
chc_m_markers<-FindAllMarkers(chc_m,only.pos = TRUE,logfc.threshold = 0.5)
chc_m_top10<-chc_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
chc_m_top10_ordergene<-chc_m_time_diff[chc_m_top10$gene,]
chc_m_time_genes<-chc_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
chc_m_time_genes<-unique(chc_m_time_genes)
chc_m_time_genes<-na.omit(chc_m_time_genes)
chc_m_newdata<-data.frame(Pseudotime = seq(min(pData(chc_m_cds)$Pseudotime),
                                           max(pData(chc_m_cds)$Pseudotime),length.out=100))
chc_m_m<-genSmoothCurves(chc_m_cds[chc_m_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = chc_m_newdata)
chc_m_m<-chc_m_m[!apply(chc_m_m,1,sum)==0,]
chc_m_m<-log10(chc_m_m+1)
chc_m_m<-chc_m_m[!apply(chc_m_m,1,sd)==0,]
chc_m_m<-Matrix::t(scale(Matrix::t(chc_m_m),center = TRUE))
chc_m_m<-chc_m_m[is.na(row.names(chc_m_m))==FALSE,]
chc_m_m[is.nan(chc_m_m)]=0
chc_m_m[chc_m_m>3]=3
chc_m_m[chc_m_m<-3]=-3
chc_m_row_dist<-as.dist((1-cor(Matrix::t(chc_m_m)))/2)
chc_m_row_dist[is.na(chc_m_row_dist)]<-1
p1<-pheatmap(chc_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = chc_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
chc_m_annotation_col<-data.frame(pseudotime=rescale(chc_m_newdata$Pseudotime,to=c(-1,1)))
row.names(chc_m_annotation_col)<-colnames(chc_m_m)
chc_m_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(chc_m_annotation_row)<-rownames(chc_m_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

chc_m_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['chc_m']]<-pheatmap(chc_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = chc_m_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = chc_m_annotation_col,
                                 annotation_colors = chc_m_anno_colors,
                                 annotation_row = chc_m_annotation_row)
main_figure[['chc_m']]<-as.ggplot(main_figure[['chc_m']])
main_figure[['chc_m']]<-ggarrange(main_figure[['chc_m']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['chc_m']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\chc_m_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

main_figure[['f3_a1']]<-DimPlot(subset(jia_t,subset = group == 'HLA-B27-'),
                                group.by = 'group',label = TRUE,cols = '#0071C2')+
  labs(title = 'HLA-B27- JIA PBMC T cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f3_a1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F3_A1.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

main_figure[['f3_a2']]<-DimPlot(subset(jia_t,subset = group == 'HLA-B27+'),
                                group.by = 'group',label = TRUE,cols = '#D75615')+
  labs(title = 'HLA-B27+ JIA PBMC T cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f3_a2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F3_A2.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

main_figure[['f3_a3']]<-DimPlot(subset(jia_t,subset = group == 'cHC'),
                                group.by = 'group',label = TRUE,cols = '#EDB11A')+
  labs(title = 'cHC PBMC T cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f3_a3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F3_A3.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

hlab27nega_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_t.rds')
hlab27nega_t_expr_matrix<-as(as.matrix(hlab27nega_t@assays$RNA@counts),'sparseMatrix')
hlab27nega_t_p_data<-hlab27nega_t@meta.data
hlab27nega_t_f_data<-data.frame(gene_short_name=row.names(hlab27nega_t),
                                row.names = row.names(hlab27nega_t))
hlab27nega_t_pd<-new('AnnotatedDataFrame',data=hlab27nega_t_p_data)
hlab27nega_t_fd<-new('AnnotatedDataFrame',data=hlab27nega_t_f_data)
hlab27nega_t_cds<-newCellDataSet(hlab27nega_t_expr_matrix,
                                 phenoData = hlab27nega_t_pd,
                                 featureData = hlab27nega_t_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27nega_t_cds<-estimateSizeFactors(hlab27nega_t_cds)
hlab27nega_t_cds<-estimateDispersions(hlab27nega_t_cds)
hlab27nega_t_cds<-detectGenes(hlab27nega_t_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27nega_t_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_t_cds<-setOrderingFilter(hlab27nega_t_cds,ordering_genes)
hlab27nega_t_cds<-reduceDimension(hlab27nega_t_cds,max_components = 2,method = 'DDRTree',num_dim = 6)
hlab27nega_t_cds<-orderCells(hlab27nega_t_cds,reverse = F)
main_figure[['f3_g1']]<-plot_cell_trajectory(hlab27nega_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  ggtitle('HLA-B27- PBMC T cell subtypes trajectory')
main_figure[['f3_g1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F3_F1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f3_g2']]<-plot_cell_trajectory(hlab27nega_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27- PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f3_g2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F3_D1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27nega_t_disp_table<-dispersionTable(hlab27nega_t_cds)
hlab27nega_t_ordering_genes<-as.character(subset(hlab27nega_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_t_time_diff<-differentialGeneTest(hlab27nega_t_cds[hlab27nega_t_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27nega_t_time_diff<-hlab27nega_t_time_diff[,c(5,2,3,4,1,6,7)]
hlab27nega_t_markers<-FindAllMarkers(hlab27nega_t,only.pos = TRUE,logfc.threshold = 0.5)
hlab27nega_t_top10<-hlab27nega_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27nega_t_top10_ordergene<-hlab27nega_t_time_diff[hlab27nega_t_top10$gene,]
hlab27nega_t_time_genes<-hlab27nega_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27nega_t_time_genes<-unique(hlab27nega_t_time_genes)
hlab27nega_t_time_genes<-na.omit(hlab27nega_t_time_genes)
hlab27nega_t_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27nega_t_cds)$Pseudotime),
                                                  max(pData(hlab27nega_t_cds)$Pseudotime),length.out=100))
hlab27nega_t_m<-genSmoothCurves(hlab27nega_t_cds[hlab27nega_t_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27nega_t_newdata)
hlab27nega_t_m<-hlab27nega_t_m[!apply(hlab27nega_t_m,1,sum)==0,]
hlab27nega_t_m<-log10(hlab27nega_t_m+1)
hlab27nega_t_m<-hlab27nega_t_m[!apply(hlab27nega_t_m,1,sd)==0,]
hlab27nega_t_m<-Matrix::t(scale(Matrix::t(hlab27nega_t_m),center = TRUE))
hlab27nega_t_m<-hlab27nega_t_m[is.na(row.names(hlab27nega_t_m))==FALSE,]
hlab27nega_t_m[is.nan(hlab27nega_t_m)]=0
hlab27nega_t_m[hlab27nega_t_m>3]=3
hlab27nega_t_m[hlab27nega_t_m<-3]=-3
hlab27nega_t_row_dist<-as.dist((1-cor(Matrix::t(hlab27nega_t_m)))/2)
hlab27nega_t_row_dist[is.na(hlab27nega_t_row_dist)]<-1
p1<-pheatmap(hlab27nega_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27nega_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
hlab27nega_t_annotation_col<-data.frame(pseudotime=rescale(hlab27nega_t_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27nega_t_annotation_col)<-colnames(hlab27nega_t_m)
hlab27nega_t_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(hlab27nega_t_annotation_row)<-rownames(hlab27nega_t_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

hlab27nega_t_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f3_h1']]<-pheatmap(hlab27nega_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27nega_t_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27nega_t_annotation_col,
                                 annotation_colors = hlab27nega_t_anno_colors,
                                 annotation_row = hlab27nega_t_annotation_row)
main_figure[['f3_h1']]<-as.ggplot(main_figure[['f3_h1']])
main_figure[['f3_h1']]<-ggarrange(main_figure[['f3_h1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f3_h1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F3_H1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

hlab27posi_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_t.rds')
hlab27posi_t_expr_matrix<-as(as.matrix(hlab27posi_t@assays$RNA@counts),'sparseMatrix')
hlab27posi_t_p_data<-hlab27posi_t@meta.data
hlab27posi_t_f_data<-data.frame(gene_short_name=row.names(hlab27posi_t),
                                row.names = row.names(hlab27posi_t))
hlab27posi_t_pd<-new('AnnotatedDataFrame',data=hlab27posi_t_p_data)
hlab27posi_t_fd<-new('AnnotatedDataFrame',data=hlab27posi_t_f_data)
hlab27posi_t_cds<-newCellDataSet(hlab27posi_t_expr_matrix,
                                 phenoData = hlab27posi_t_pd,
                                 featureData = hlab27posi_t_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27posi_t_cds<-estimateSizeFactors(hlab27posi_t_cds)
hlab27posi_t_cds<-estimateDispersions(hlab27posi_t_cds)
hlab27posi_t_cds<-detectGenes(hlab27posi_t_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27posi_t_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_t_cds<-setOrderingFilter(hlab27posi_t_cds,ordering_genes)
hlab27posi_t_cds<-reduceDimension(hlab27posi_t_cds,max_components = 2,method = 'DDRTree',num_dim = 6)
hlab27posi_t_cds<-orderCells(hlab27posi_t_cds,reverse = F)
main_figure[['f3_g3']]<-plot_cell_trajectory(hlab27posi_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#01C4C7','#55BBE1','#90E2F3','#CFF0F7'))+
  ggtitle('HLA-B27+ PBMC T cell subtypes trajectory')
main_figure[['f3_g3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F3_D1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f3_g4']]<-plot_cell_trajectory(hlab27posi_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27+ PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f3_g4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F3_G4.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27posi_t_disp_table<-dispersionTable(hlab27posi_t_cds)
hlab27posi_t_ordering_genes<-as.character(subset(hlab27posi_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_t_time_diff<-differentialGeneTest(hlab27posi_t_cds[hlab27posi_t_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27posi_t_time_diff<-hlab27posi_t_time_diff[,c(5,2,3,4,1,6,7)]
hlab27posi_t_markers<-FindAllMarkers(hlab27posi_t,only.pos = TRUE,logfc.threshold = 0.5)
hlab27posi_t_top10<-hlab27posi_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27posi_t_top10_ordergene<-hlab27posi_t_time_diff[hlab27posi_t_top10$gene,]
hlab27posi_t_time_genes<-hlab27posi_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27posi_t_time_genes<-unique(hlab27posi_t_time_genes)
hlab27posi_t_time_genes<-na.omit(hlab27posi_t_time_genes)
hlab27posi_t_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27posi_t_cds)$Pseudotime),
                                                  max(pData(hlab27posi_t_cds)$Pseudotime),length.out=100))
hlab27posi_t_m<-genSmoothCurves(hlab27posi_t_cds[hlab27posi_t_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27posi_t_newdata)
hlab27posi_t_m<-hlab27posi_t_m[!apply(hlab27posi_t_m,1,sum)==0,]
hlab27posi_t_m<-log10(hlab27posi_t_m+1)
hlab27posi_t_m<-hlab27posi_t_m[!apply(hlab27posi_t_m,1,sd)==0,]
hlab27posi_t_m<-Matrix::t(scale(Matrix::t(hlab27posi_t_m),center = TRUE))
hlab27posi_t_m<-hlab27posi_t_m[is.na(row.names(hlab27posi_t_m))==FALSE,]
hlab27posi_t_m[is.nan(hlab27posi_t_m)]=0
hlab27posi_t_m[hlab27posi_t_m>3]=3
hlab27posi_t_m[hlab27posi_t_m<-3]=-3
hlab27posi_t_row_dist<-as.dist((1-cor(Matrix::t(hlab27posi_t_m)))/2)
hlab27posi_t_row_dist[is.na(hlab27posi_t_row_dist)]<-1
p1<-pheatmap(hlab27posi_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27posi_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
hlab27posi_t_annotation_col<-data.frame(pseudotime=rescale(hlab27posi_t_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27posi_t_annotation_col)<-colnames(hlab27posi_t_m)
hlab27posi_t_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(hlab27posi_t_annotation_row)<-rownames(hlab27posi_t_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

hlab27posi_t_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f3_h2']]<-pheatmap(hlab27posi_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27posi_t_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27posi_t_annotation_col,
                                 annotation_colors = hlab27posi_t_anno_colors,
                                 annotation_row = hlab27posi_t_annotation_row)
main_figure[['f3_h2']]<-as.ggplot(main_figure[['f3_h2']])
main_figure[['f3_h2']]<-ggarrange(main_figure[['f3_h2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f3_h2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F3_H2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

main_figure[['f4_a1']]<-DimPlot(subset(jia_b,subset = group == 'HLA-B27-'),
                                group.by = 'group',label = TRUE,cols = '#0071C2')+
  labs(title = 'HLA-B27- JIA PBMC B cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f4_a1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F4_A1.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

main_figure[['f4_a2']]<-DimPlot(subset(jia_b,subset = group == 'HLA-B27+'),
                                group.by = 'group',label = TRUE,cols = '#D75615')+
  labs(title = 'HLA-B27+ JIA PBMC B cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f4_a2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F4_A2.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

main_figure[['f4_a3']]<-DimPlot(subset(jia_b,subset = group == 'cHC'),
                                group.by = 'group',label = TRUE,cols = '#EDB11A')+
  labs(title = 'cHC PBMC B cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f4_a3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F4_A3.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

hlab27nega_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_b.rds')
hlab27nega_b_expr_matrix<-as(as.matrix(hlab27nega_b@assays$RNA@counts),'sparseMatrix')
hlab27nega_b_p_data<-hlab27nega_b@meta.data
hlab27nega_b_f_data<-data.frame(gene_short_name=row.names(hlab27nega_b),
                                row.names = row.names(hlab27nega_b))
hlab27nega_b_pd<-new('AnnotatedDataFrame',data=hlab27nega_b_p_data)
hlab27nega_b_fd<-new('AnnotatedDataFrame',data=hlab27nega_b_f_data)
hlab27nega_b_cds<-newCellDataSet(hlab27nega_b_expr_matrix,
                                 phenoData = hlab27nega_b_pd,
                                 featureData = hlab27nega_b_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27nega_b_cds<-estimateSizeFactors(hlab27nega_b_cds)
hlab27nega_b_cds<-estimateDispersions(hlab27nega_b_cds)
hlab27nega_b_cds<-detectGenes(hlab27nega_b_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27nega_b_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_b_cds<-setOrderingFilter(hlab27nega_b_cds,ordering_genes)
hlab27nega_b_cds<-reduceDimension(hlab27nega_b_cds,max_components = 2,method = 'DDRTree')
hlab27nega_b_cds<-orderCells(hlab27nega_b_cds,reverse = F)
main_figure[['f4_f1']]<-plot_cell_trajectory(hlab27nega_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#EC8D63','#DE4B3F','#DE4247'))+
  ggtitle('HLA-B27- PBMC B cell subtypes trajectory')
main_figure[['f4_f1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_O1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f4_f2']]<-plot_cell_trajectory(hlab27nega_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27- PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f4_f2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_O2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27nega_b_disp_table<-dispersionTable(hlab27nega_b_cds)
hlab27nega_b_ordering_genes<-as.character(subset(hlab27nega_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_b_time_diff<-differentialGeneTest(hlab27nega_b_cds[hlab27nega_b_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27nega_b_time_diff<-hlab27nega_b_time_diff[,c(5,2,3,4,1,6,7)]
hlab27nega_b_markers<-FindAllMarkers(hlab27nega_b,only.pos = TRUE,logfc.threshold = 0.5)
hlab27nega_b_top10<-hlab27nega_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27nega_b_top10_ordergene<-hlab27nega_b_time_diff[hlab27nega_b_top10$gene,]
hlab27nega_b_time_genes<-hlab27nega_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27nega_b_time_genes<-unique(hlab27nega_b_time_genes)
hlab27nega_b_time_genes<-na.omit(hlab27nega_b_time_genes)
hlab27nega_b_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27nega_b_cds)$Pseudotime),
                                                  max(pData(hlab27nega_b_cds)$Pseudotime),length.out=100))
hlab27nega_b_m<-genSmoothCurves(hlab27nega_b_cds[hlab27nega_b_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27nega_b_newdata)
hlab27nega_b_m<-hlab27nega_b_m[!apply(hlab27nega_b_m,1,sum)==0,]
hlab27nega_b_m<-log10(hlab27nega_b_m+1)
hlab27nega_b_m<-hlab27nega_b_m[!apply(hlab27nega_b_m,1,sd)==0,]
hlab27nega_b_m<-Matrix::t(scale(Matrix::t(hlab27nega_b_m),center = TRUE))
hlab27nega_b_m<-hlab27nega_b_m[is.na(row.names(hlab27nega_b_m))==FALSE,]
hlab27nega_b_m[is.nan(hlab27nega_b_m)]=0
hlab27nega_b_m[hlab27nega_b_m>3]=3
hlab27nega_b_m[hlab27nega_b_m<-3]=-3
hlab27nega_b_row_dist<-as.dist((1-cor(Matrix::t(hlab27nega_b_m)))/2)
hlab27nega_b_row_dist[is.na(hlab27nega_b_row_dist)]<-1

p1<-pheatmap(hlab27nega_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27nega_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1

hlab27nega_b_annotation_col<-data.frame(pseudotime=rescale(hlab27nega_b_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27nega_b_annotation_col)<-colnames(hlab27nega_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(hlab27nega_b_m)
hlab27nega_b_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f4_g1']]<-pheatmap(hlab27nega_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27nega_b_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27nega_b_annotation_col,
                                 annotation_colors = hlab27nega_b_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f4_g1']]<-as.ggplot(main_figure[['f4_g1']])
main_figure[['f4_g1']]<-ggarrange(main_figure[['f4_g1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f4_g1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F4_G1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

hlab27posi_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_b.rds')
hlab27posi_b_expr_matrix<-as(as.matrix(hlab27posi_b@assays$RNA@counts),'sparseMatrix')
hlab27posi_b_p_data<-hlab27posi_b@meta.data
hlab27posi_b_f_data<-data.frame(gene_short_name=row.names(hlab27posi_b),
                                row.names = row.names(hlab27posi_b))
hlab27posi_b_pd<-new('AnnotatedDataFrame',data=hlab27posi_b_p_data)
hlab27posi_b_fd<-new('AnnotatedDataFrame',data=hlab27posi_b_f_data)
hlab27posi_b_cds<-newCellDataSet(hlab27posi_b_expr_matrix,
                                 phenoData = hlab27posi_b_pd,
                                 featureData = hlab27posi_b_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27posi_b_cds<-estimateSizeFactors(hlab27posi_b_cds)
hlab27posi_b_cds<-estimateDispersions(hlab27posi_b_cds)
hlab27posi_b_cds<-detectGenes(hlab27posi_b_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27posi_b_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_b_cds<-setOrderingFilter(hlab27posi_b_cds,ordering_genes)
hlab27posi_b_cds<-reduceDimension(hlab27posi_b_cds,max_components = 2,method = 'DDRTree')
hlab27posi_b_cds<-orderCells(hlab27posi_b_cds,reverse = F)
main_figure[['f4_f3']]<-plot_cell_trajectory(hlab27posi_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#EC8D63','#DE4B3F','#DE4247'))+
  ggtitle('HLA-B27+ PBMC B cell subtypes trajectory')
main_figure[['f4_f3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F5_M1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f4_f4']]<-plot_cell_trajectory(hlab27posi_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27+ PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f4_f4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F4_F4.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27posi_b_disp_table<-dispersionTable(hlab27posi_b_cds)
hlab27posi_b_ordering_genes<-as.character(subset(hlab27posi_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_b_time_diff<-differentialGeneTest(hlab27posi_b_cds[hlab27posi_b_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27posi_b_time_diff<-hlab27posi_b_time_diff[,c(5,2,3,4,1,6,7)]
hlab27posi_b_markers<-FindAllMarkers(hlab27posi_b,only.pos = TRUE,logfc.threshold = 0.5)
hlab27posi_b_top10<-hlab27posi_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27posi_b_top10_ordergene<-hlab27posi_b_time_diff[hlab27posi_b_top10$gene,]
hlab27posi_b_time_genes<-hlab27posi_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27posi_b_time_genes<-unique(hlab27posi_b_time_genes)
hlab27posi_b_time_genes<-na.omit(hlab27posi_b_time_genes)
hlab27posi_b_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27posi_b_cds)$Pseudotime),
                                                  max(pData(hlab27posi_b_cds)$Pseudotime),length.out=100))
hlab27posi_b_m<-genSmoothCurves(hlab27posi_b_cds[hlab27posi_b_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27posi_b_newdata)
hlab27posi_b_m<-hlab27posi_b_m[!apply(hlab27posi_b_m,1,sum)==0,]
hlab27posi_b_m<-log10(hlab27posi_b_m+1)
hlab27posi_b_m<-hlab27posi_b_m[!apply(hlab27posi_b_m,1,sd)==0,]
hlab27posi_b_m<-Matrix::t(scale(Matrix::t(hlab27posi_b_m),center = TRUE))
hlab27posi_b_m<-hlab27posi_b_m[is.na(row.names(hlab27posi_b_m))==FALSE,]
hlab27posi_b_m[is.nan(hlab27posi_b_m)]=0
hlab27posi_b_m[hlab27posi_b_m>3]=3
hlab27posi_b_m[hlab27posi_b_m<-3]=-3
hlab27posi_b_row_dist<-as.dist((1-cor(Matrix::t(hlab27posi_b_m)))/2)
hlab27posi_b_row_dist[is.na(hlab27posi_b_row_dist)]<-1

p1<-pheatmap(hlab27posi_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27posi_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1

hlab27posi_b_annotation_col<-data.frame(pseudotime=rescale(hlab27posi_b_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27posi_b_annotation_col)<-colnames(hlab27posi_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(hlab27posi_b_m)
hlab27posi_b_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f4_g2']]<-pheatmap(hlab27posi_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27posi_b_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27posi_b_annotation_col,
                                 annotation_colors = hlab27posi_b_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f4_g2']]<-as.ggplot(main_figure[['f4_g2']])
main_figure[['f4_g2']]<-ggarrange(main_figure[['f4_g2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f4_g2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F4_G2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

main_figure[['f5_a1']]<-DimPlot(subset(jia_m,subset = group == 'HLA-B27-'),
                                group.by = 'group',label = TRUE,cols = '#0071C2')+
  labs(title = 'HLA-B27- JIA PBMC Myeloid cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f5_a1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F5_A1.pdf',units = 'in',dpi = 300,
       width = 7,height = 6,device = "pdf")

main_figure[['f5_a2']]<-DimPlot(subset(jia_m,subset = group == 'HLA-B27+'),
                                group.by = 'group',label = TRUE,cols = '#D75615')+
  labs(title = 'HLA-B27+ JIA PBMC Myeloid cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f5_a2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F5_A2.pdf',units = 'in',dpi = 300,
       width = 7,height = 6,device = "pdf")

main_figure[['f5_a3']]<-DimPlot(subset(jia_m,subset = group == 'cHC'),
                                group.by = 'group',label = TRUE,cols = '#EDB11A')+
  labs(title = 'cHC PBMC Myeloid cells distribution',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f5_a3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure\\F5_A3.pdf',units = 'in',dpi = 300,
       width = 7,height = 6,device = "pdf")

hlab27nega_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27-\\hlab27-_m.rds')
hlab27nega_m<-subset(hlab27nega_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
hlab27nega_m_expr_matrix<-as(as.matrix(hlab27nega_m@assays$RNA@counts),'sparseMatrix')
hlab27nega_m_p_data<-hlab27nega_m@meta.data
hlab27nega_m_f_data<-data.frame(gene_short_name=row.names(hlab27nega_m),
                                row.names = row.names(hlab27nega_m))
hlab27nega_m_pd<-new('AnnotatedDataFrame',data=hlab27nega_m_p_data)
hlab27nega_m_fd<-new('AnnotatedDataFrame',data=hlab27nega_m_f_data)
hlab27nega_m_cds<-newCellDataSet(hlab27nega_m_expr_matrix,
                                 phenoData = hlab27nega_m_pd,
                                 featureData = hlab27nega_m_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27nega_m_cds<-estimateSizeFactors(hlab27nega_m_cds)
hlab27nega_m_cds<-estimateDispersions(hlab27nega_m_cds)
hlab27nega_m_cds<-detectGenes(hlab27nega_m_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27nega_m_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_m_cds<-setOrderingFilter(hlab27nega_m_cds,ordering_genes)
hlab27nega_m_cds<-reduceDimension(hlab27nega_m_cds,max_components = 2,method = 'DDRTree')
hlab27nega_m_cds<-orderCells(hlab27nega_m_cds,reverse = F)
main_figure[['f5_f1']]<-plot_cell_trajectory(hlab27nega_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#65A644','#657A51','#8DC591'))+
  ggtitle('HLA-B27- PBMC Monocyte subtypes trajectory')
main_figure[['f5_f1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F6_O1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f5_f2']]<-plot_cell_trajectory(hlab27nega_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27- PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f5_f2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F5_F2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27nega_m_disp_table<-dispersionTable(hlab27nega_m_cds)
hlab27nega_m_ordering_genes<-as.character(subset(hlab27nega_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27nega_m_time_diff<-differentialGeneTest(hlab27nega_m_cds[hlab27nega_m_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27nega_m_time_diff<-hlab27nega_m_time_diff[,c(5,2,3,4,1,6,7)]
hlab27nega_m_markers<-FindAllMarkers(hlab27nega_m,only.pos = TRUE,logfc.threshold = 0.5)
hlab27nega_m_top10<-hlab27nega_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27nega_m_top10_ordergene<-hlab27nega_m_time_diff[hlab27nega_m_top10$gene,]
hlab27nega_m_time_genes<-hlab27nega_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27nega_m_time_genes<-unique(hlab27nega_m_time_genes)
hlab27nega_m_time_genes<-na.omit(hlab27nega_m_time_genes)
hlab27nega_m_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27nega_m_cds)$Pseudotime),
                                                  max(pData(hlab27nega_m_cds)$Pseudotime),length.out=100))
hlab27nega_m_m<-genSmoothCurves(hlab27nega_m_cds[hlab27nega_m_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27nega_m_newdata)
hlab27nega_m_m<-hlab27nega_m_m[!apply(hlab27nega_m_m,1,sum)==0,]
hlab27nega_m_m<-log10(hlab27nega_m_m+1)
hlab27nega_m_m<-hlab27nega_m_m[!apply(hlab27nega_m_m,1,sd)==0,]
hlab27nega_m_m<-Matrix::t(scale(Matrix::t(hlab27nega_m_m),center = TRUE))
hlab27nega_m_m<-hlab27nega_m_m[is.na(row.names(hlab27nega_m_m))==FALSE,]
hlab27nega_m_m[is.nan(hlab27nega_m_m)]=0
hlab27nega_m_m[hlab27nega_m_m>3]=3
hlab27nega_m_m[hlab27nega_m_m<-3]=-3
hlab27nega_m_row_dist<-as.dist((1-cor(Matrix::t(hlab27nega_m_m)))/2)
hlab27nega_m_row_dist[is.na(hlab27nega_m_row_dist)]<-1

p1<-pheatmap(hlab27nega_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27nega_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1

hlab27nega_m_annotation_col<-data.frame(pseudotime=rescale(hlab27nega_m_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27nega_m_annotation_col)<-colnames(hlab27nega_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(hlab27nega_m_m)
hlab27nega_m_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f5_g1']]<-pheatmap(hlab27nega_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27nega_m_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27nega_m_annotation_col,
                                 annotation_colors = hlab27nega_m_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f5_g1']]<-as.ggplot(main_figure[['f5_g1']])
main_figure[['f5_g1']]<-ggarrange(main_figure[['f5_g1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f5_g1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F5_G1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

hlab27posi_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\HLA-B27+\\hlab27+_m.rds')
hlab27posi_m<-subset(hlab27posi_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
hlab27posi_m_expr_matrix<-as(as.matrix(hlab27posi_m@assays$RNA@counts),'sparseMatrix')
hlab27posi_m_p_data<-hlab27posi_m@meta.data
hlab27posi_m_f_data<-data.frame(gene_short_name=row.names(hlab27posi_m),
                                row.names = row.names(hlab27posi_m))
hlab27posi_m_pd<-new('AnnotatedDataFrame',data=hlab27posi_m_p_data)
hlab27posi_m_fd<-new('AnnotatedDataFrame',data=hlab27posi_m_f_data)
hlab27posi_m_cds<-newCellDataSet(hlab27posi_m_expr_matrix,
                                 phenoData = hlab27posi_m_pd,
                                 featureData = hlab27posi_m_fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
hlab27posi_m_cds<-estimateSizeFactors(hlab27posi_m_cds)
hlab27posi_m_cds<-estimateDispersions(hlab27posi_m_cds)
hlab27posi_m_cds<-detectGenes(hlab27posi_m_cds,min_expr = 0.5)
disp_table<-dispersionTable(hlab27posi_m_cds)
ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_m_cds<-setOrderingFilter(hlab27posi_m_cds,ordering_genes)
hlab27posi_m_cds<-reduceDimension(hlab27posi_m_cds,max_components = 2,method = 'DDRTree')
hlab27posi_m_cds<-orderCells(hlab27posi_m_cds,reverse = F)
main_figure[['f5_f3']]<-plot_cell_trajectory(hlab27posi_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#65A644','#657A51','#8DC591'))+
  ggtitle('HLA-B27+ PBMC Monocyte subtypes trajectory')
main_figure[['f5_f3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F6_M1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f5_f4']]<-plot_cell_trajectory(hlab27posi_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('HLA-B27+ PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f5_f4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F5_F4.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

hlab27posi_m_disp_table<-dispersionTable(hlab27posi_m_cds)
hlab27posi_m_ordering_genes<-as.character(subset(hlab27posi_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
hlab27posi_m_time_diff<-differentialGeneTest(hlab27posi_m_cds[hlab27posi_m_ordering_genes,],cores = 1,
                                             fullModelFormulaStr = '~sm.ns(Pseudotime)')
hlab27posi_m_time_diff<-hlab27posi_m_time_diff[,c(5,2,3,4,1,6,7)]
hlab27posi_m_markers<-FindAllMarkers(hlab27posi_m,only.pos = TRUE,logfc.threshold = 0.5)
hlab27posi_m_top10<-hlab27posi_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
hlab27posi_m_top10_ordergene<-hlab27posi_m_time_diff[hlab27posi_m_top10$gene,]
hlab27posi_m_time_genes<-hlab27posi_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
hlab27posi_m_time_genes<-unique(hlab27posi_m_time_genes)
hlab27posi_m_time_genes<-na.omit(hlab27posi_m_time_genes)
hlab27posi_m_newdata<-data.frame(Pseudotime = seq(min(pData(hlab27posi_m_cds)$Pseudotime),
                                                  max(pData(hlab27posi_m_cds)$Pseudotime),length.out=100))
hlab27posi_m_m<-genSmoothCurves(hlab27posi_m_cds[hlab27posi_m_time_genes],
                                trend_formula = '~sm.ns(Pseudotime,df=3)',
                                relative_expr = T,new_data = hlab27posi_m_newdata)
hlab27posi_m_m<-hlab27posi_m_m[!apply(hlab27posi_m_m,1,sum)==0,]
hlab27posi_m_m<-log10(hlab27posi_m_m+1)
hlab27posi_m_m<-hlab27posi_m_m[!apply(hlab27posi_m_m,1,sd)==0,]
hlab27posi_m_m<-Matrix::t(scale(Matrix::t(hlab27posi_m_m),center = TRUE))
hlab27posi_m_m<-hlab27posi_m_m[is.na(row.names(hlab27posi_m_m))==FALSE,]
hlab27posi_m_m[is.nan(hlab27posi_m_m)]=0
hlab27posi_m_m[hlab27posi_m_m>3]=3
hlab27posi_m_m[hlab27posi_m_m<-3]=-3
hlab27posi_m_row_dist<-as.dist((1-cor(Matrix::t(hlab27posi_m_m)))/2)
hlab27posi_m_row_dist[is.na(hlab27posi_m_row_dist)]<-1

p1<-pheatmap(hlab27posi_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = hlab27posi_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1

hlab27posi_m_annotation_col<-data.frame(pseudotime=rescale(hlab27posi_m_newdata$Pseudotime,to=c(-1,1)))
row.names(hlab27posi_m_annotation_col)<-colnames(hlab27posi_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(hlab27posi_m_m)
hlab27posi_m_anno_colors<-list(pseudotime=viridis(100),
                               Cluster=row_color)
main_figure[['f5_g2']]<-pheatmap(hlab27posi_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = hlab27posi_m_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = hlab27posi_m_annotation_col,
                                 annotation_colors = hlab27posi_m_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f5_g2']]<-as.ggplot(main_figure[['f5_g2']])
main_figure[['f5_g2']]<-ggarrange(main_figure[['f5_g2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f5_g2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F5_G2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

pss_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\pSS\\pss_t.rds')
pss_t<-NormalizeData(pss_t)
pss_t<-FindVariableFeatures(pss_t,selection.method = 'vst',nfeatures = 2000)
pss_t<-ScaleData(pss_t)
pss_t<-RunPCA(pss_t,features = VariableFeatures(object = pss_t))
ElbowPlot(pss_t)
pss_t<-FindNeighbors(pss_t,dims = 1:15)
pss_t<-FindClusters(pss_t,resolution = 0.3)
pss_t<-RunUMAP(pss_t,dims = 1:15)
main_figure[['f6_g1']]<-DimPlot(pss_t,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(7,'Set1')))+
  labs(title = 'pSS PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_g1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_G1.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

pss_t_expr_matrix<-as(as.matrix(pss_t@assays$RNA@counts),'sparseMatrix')
pss_t_p_data<-pss_t@meta.data
pss_t_f_data<-data.frame(gene_short_name=row.names(pss_t),row.names = row.names(pss_t))
pss_t_pd<-new('AnnotatedDataFrame',data=pss_t_p_data)
pss_t_fd<-new('AnnotatedDataFrame',data=pss_t_f_data)
pss_t_cds<-newCellDataSet(pss_t_expr_matrix,
                          phenoData = pss_t_pd,
                          featureData = pss_t_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
pss_t_cds<-estimateSizeFactors(pss_t_cds)
pss_t_cds<-estimateDispersions(pss_t_cds)
pss_t_cds<-detectGenes(pss_t_cds,min_expr = 0.5)
pss_t_disp_table<-dispersionTable(pss_t_cds)
pss_t_ordering_genes<-as.character(subset(pss_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
pss_t_cds<-setOrderingFilter(pss_t_cds,pss_t_ordering_genes)
pss_t_cds<-reduceDimension(pss_t_cds,max_components = 2,method = 'DDRTree')
pss_t_cds<-orderCells(pss_t_cds,reverse = F)
main_figure[['f6_h1']]<-plot_cell_trajectory(pss_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#0047AB','#01C4C7','#55BBE1','#90E2F3','#007BA7'))+
  ggtitle('pSS PBMC T cell subtypes trajectory')
main_figure[['f6_h1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F4_H1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_h2']]<-plot_cell_trajectory(pss_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('pSS PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f6_h2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_H2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

#pss_t_disp_table<-dispersionTable(pss_t_cds)
#jia_t_ordering_genes<-as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
pss_t_time_diff<-differentialGeneTest(pss_t_cds[pss_t_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
pss_t_time_diff<-pss_t_time_diff[,c(5,2,3,4,1,6,7)]
pss_t_markers<-FindAllMarkers(pss_t,only.pos = TRUE,logfc.threshold = 0.5)
pss_t_top10<-pss_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
pss_t_top10_ordergene<-pss_t_time_diff[pss_t_top10$gene,]
pss_t_time_genes<-pss_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
pss_t_time_genes<-unique(pss_t_time_genes)
pss_t_time_genes<-na.omit(pss_t_time_genes)
pss_t_newdata<-data.frame(Pseudotime = seq(min(pData(pss_t_cds)$Pseudotime),
                                           max(pData(pss_t_cds)$Pseudotime),length.out=100))
pss_t_m<-genSmoothCurves(pss_t_cds[pss_t_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = pss_t_newdata)
pss_t_m<-pss_t_m[!apply(pss_t_m,1,sum)==0,]
pss_t_m<-log10(pss_t_m+1)
pss_t_m<-pss_t_m[!apply(pss_t_m,1,sd)==0,]
pss_t_m<-Matrix::t(scale(Matrix::t(pss_t_m),center = TRUE))
pss_t_m<-pss_t_m[is.na(row.names(pss_t_m))==FALSE,]
pss_t_m[is.nan(pss_t_m)]=0
pss_t_m[pss_t_m>3]=3
pss_t_m[pss_t_m<-3]=-3
pss_t_row_dist<-as.dist((1-cor(Matrix::t(pss_t_m)))/2)
pss_t_row_dist[is.na(pss_t_row_dist)]<-1

p1<-pheatmap(pss_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = pss_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

pss_t_annotation_col<-data.frame(pseudotime=rescale(pss_t_newdata$Pseudotime,to=c(-1,1)))
row.names(pss_t_annotation_col)<-colnames(pss_t_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(pss_t_m)
pss_t_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_i1']]<-pheatmap(pss_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = pss_t_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = pss_t_annotation_col,
                                 annotation_colors = pss_t_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_i1']]<-as.ggplot(main_figure[['f6_i1']])
main_figure[['f6_i1']]<-ggarrange(main_figure[['f6_i1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_i1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_I1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

sle_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\SLE\\sle_t.rds')
sle_t<-NormalizeData(sle_t)
sle_t<-FindVariableFeatures(sle_t,selection.method = 'vst',nfeatures = 2000)
sle_t<-ScaleData(sle_t)
sle_t<-RunPCA(sle_t,features = VariableFeatures(object = sle_t))
ElbowPlot(sle_t)
sle_t<-FindNeighbors(sle_t,dims = 1:15)
sle_t<-FindClusters(sle_t,resolution = 0.3)
sle_t<-RunUMAP(sle_t,dims = 1:15)
main_figure[['f6_g2']]<-DimPlot(sle_t,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(5,'Set1')))+
  labs(title = 'SLE PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_g2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_G2.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

sle_t_expr_matrix<-as(as.matrix(sle_t@assays$RNA@counts),'sparseMatrix')
sle_t_p_data<-sle_t@meta.data
sle_t_f_data<-data.frame(gene_short_name=row.names(sle_t),row.names = row.names(sle_t))
sle_t_pd<-new('AnnotatedDataFrame',data=sle_t_p_data)
sle_t_fd<-new('AnnotatedDataFrame',data=sle_t_f_data)
sle_t_cds<-newCellDataSet(sle_t_expr_matrix,
                          phenoData = sle_t_pd,
                          featureData = sle_t_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
sle_t_cds<-estimateSizeFactors(sle_t_cds)
sle_t_cds<-estimateDispersions(sle_t_cds)
sle_t_cds<-detectGenes(sle_t_cds,min_expr = 0.5)
sle_t_disp_table<-dispersionTable(sle_t_cds)
sle_t_ordering_genes<-as.character(subset(sle_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
sle_t_cds<-setOrderingFilter(sle_t_cds,sle_t_ordering_genes)
sle_t_cds<-reduceDimension(sle_t_cds,max_components = 2,method = 'DDRTree')
sle_t_cds<-orderCells(sle_t_cds,reverse = F)
main_figure[['f6_h3']]<-plot_cell_trajectory(sle_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  scale_color_manual(values = c('#248CBF','#3C9FCD','#0047AB','#01C4C7','#55BBE1'))+
  ggtitle('SLE PBMC T cell subtypes trajectory')
main_figure[['f6_h3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F4_J1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_h4']]<-plot_cell_trajectory(sle_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('SLE PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f6_h4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-9-12\\F4_J2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

sle_t_time_diff<-differentialGeneTest(sle_t_cds[sle_t_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
sle_t_time_diff<-sle_t_time_diff[,c(5,2,3,4,1,6,7)]
sle_t_markers<-FindAllMarkers(sle_t,only.pos = TRUE,logfc.threshold = 0.5)
sle_t_top10<-sle_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
sle_t_top10_ordergene<-sle_t_time_diff[sle_t_top10$gene,]
sle_t_time_genes<-sle_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
sle_t_time_genes<-unique(sle_t_time_genes)
sle_t_time_genes<-na.omit(sle_t_time_genes)
sle_t_newdata<-data.frame(Pseudotime = seq(min(pData(sle_t_cds)$Pseudotime),
                                           max(pData(sle_t_cds)$Pseudotime),length.out=100))
sle_t_m<-genSmoothCurves(sle_t_cds[sle_t_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = sle_t_newdata)
sle_t_m<-sle_t_m[!apply(sle_t_m,1,sum)==0,]
sle_t_m<-log10(sle_t_m+1)
sle_t_m<-sle_t_m[!apply(sle_t_m,1,sd)==0,]
sle_t_m<-Matrix::t(scale(Matrix::t(sle_t_m),center = TRUE))
sle_t_m<-sle_t_m[is.na(row.names(sle_t_m))==FALSE,]
sle_t_m[is.nan(sle_t_m)]=0
sle_t_m[sle_t_m>3]=3
sle_t_m[sle_t_m<-3]=-3
sle_t_row_dist<-as.dist((1-cor(Matrix::t(sle_t_m)))/2)
sle_t_row_dist[is.na(sle_t_row_dist)]<-1

p1<-pheatmap(sle_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = sle_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

sle_t_annotation_col<-data.frame(pseudotime=rescale(sle_t_newdata$Pseudotime,to=c(-1,1)))
row.names(sle_t_annotation_col)<-colnames(sle_t_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(sle_t_m)
sle_t_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_i2']]<-pheatmap(sle_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,clustering_method = 'ward.D2',
                                 clustering_distance_rows = sle_t_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = sle_t_annotation_col,
                                 annotation_colors = sle_t_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_i2']]<-as.ggplot(main_figure[['f6_i2']])
main_figure[['f6_i2']]<-ggarrange(main_figure[['f6_i2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_i2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_I2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

ahc_t<-readRDS('D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_t.rds')
ahc_t<-NormalizeData(ahc_t)
ahc_t<-FindVariableFeatures(ahc_t,selection.method = 'vst',nfeatures = 2000)
ahc_t<-ScaleData(ahc_t)
ahc_t<-RunPCA(ahc_t,features = VariableFeatures(object = ahc_t))
ElbowPlot(ahc_t)
ahc_t<-FindNeighbors(ahc_t,dims = 1:15)
ahc_t<-FindClusters(ahc_t,resolution = 0.3)
ahc_t<-RunUMAP(ahc_t,dims = 1:15)
main_figure[['f6_g3']]<-DimPlot(ahc_t,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(8,'Set1')))+
  labs(title = 'aHC PBMC T cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_g3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_G3.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

ahc_t_expr_matrix<-as(as.matrix(ahc_t@assays$RNA@counts),'sparseMatrix')
ahc_t_p_data<-ahc_t@meta.data
ahc_t_f_data<-data.frame(gene_short_name=row.names(ahc_t),row.names = row.names(ahc_t))
ahc_t_pd<-new('AnnotatedDataFrame',data=ahc_t_p_data)
ahc_t_fd<-new('AnnotatedDataFrame',data=ahc_t_f_data)
ahc_t_cds<-newCellDataSet(ahc_t_expr_matrix,
                          phenoData = ahc_t_pd,
                          featureData = ahc_t_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
ahc_t_cds<-estimateSizeFactors(ahc_t_cds)
ahc_t_cds<-estimateDispersions(ahc_t_cds)
ahc_t_cds<-detectGenes(ahc_t_cds,min_expr = 0.5)
ahc_t_disp_table<-dispersionTable(ahc_t_cds)
ahc_t_ordering_genes<-as.character(subset(ahc_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
ahc_t_cds<-setOrderingFilter(ahc_t_cds,ahc_t_ordering_genes)
ahc_t_cds<-reduceDimension(ahc_t_cds,max_components = 2,method = 'DDRTree')
ahc_t_cds<-orderCells(ahc_t_cds,reverse = F)
main_figure[['f6_h5']]<-plot_cell_trajectory(ahc_t_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('aHC PBMC T cell subtypes trajectory')
main_figure[['f6_h5']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_H5.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_h6']]<-plot_cell_trajectory(ahc_t_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('aHC PBMC T cell subtypes trajectory (Pseudotime)')
main_figure[['f6_h6']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_H6.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

ahc_t_disp_table<-dispersionTable(ahc_t_cds)
ahc_t_ordering_genes<-as.character(subset(ahc_t_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
ahc_t_time_diff<-differentialGeneTest(ahc_t_cds[ahc_t_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
ahc_t_time_diff<-ahc_t_time_diff[,c(5,2,3,4,1,6,7)]
ahc_t_markers<-FindAllMarkers(ahc_t,only.pos = TRUE,logfc.threshold = 0.5)
ahc_t_top10<-ahc_t_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
ahc_t_top10_ordergene<-ahc_t_time_diff[ahc_t_top10$gene,]
ahc_t_time_genes<-ahc_t_top10_ordergene %>% pull(gene_short_name) %>% as.character()
ahc_t_time_genes<-unique(ahc_t_time_genes)
ahc_t_time_genes<-na.omit(ahc_t_time_genes)
ahc_t_newdata<-data.frame(Pseudotime = seq(min(pData(ahc_t_cds)$Pseudotime),
                                           max(pData(ahc_t_cds)$Pseudotime),length.out=100))
ahc_t_m<-genSmoothCurves(ahc_t_cds[ahc_t_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = ahc_t_newdata)
ahc_t_m<-ahc_t_m[!apply(ahc_t_m,1,sum)==0,]
ahc_t_m<-log10(ahc_t_m+1)
ahc_t_m<-ahc_t_m[!apply(ahc_t_m,1,sd)==0,]
ahc_t_m<-Matrix::t(scale(Matrix::t(ahc_t_m),center = TRUE))
ahc_t_m<-ahc_t_m[is.na(row.names(ahc_t_m))==FALSE,]
ahc_t_m[is.nan(ahc_t_m)]=0
ahc_t_m[ahc_t_m>3]=3
ahc_t_m[ahc_t_m<-3]=-3
ahc_t_row_dist<-as.dist((1-cor(Matrix::t(ahc_t_m)))/2)
ahc_t_row_dist[is.na(ahc_t_row_dist)]<-1
p1<-pheatmap(ahc_t_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = ahc_t_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))
p1
ahc_t_annotation_col<-data.frame(pseudotime=rescale(ahc_t_newdata$Pseudotime,to=c(-1,1)))
row.names(ahc_t_annotation_col)<-colnames(ahc_t_m)
ahc_t_annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(ahc_t_annotation_row)<-rownames(ahc_t_m)
row_color<-c('#85B22E','#E29827','#922927','#57C3F3')
names(row_color)<-c('1','2','3','4')

ahc_t_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['ahc_t']]<-pheatmap(ahc_t_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = ahc_t_row_dist,
                                 border_color = NA,filename = NA,clustering_method = 'ward.D2',
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = ahc_t_annotation_col,
                                 annotation_colors = ahc_t_anno_colors,
                                 annotation_row = ahc_t_annotation_row)
main_figure[['ahc_t']]<-as.ggplot(main_figure[['ahc_t']])
main_figure[['ahc_t']]<-ggarrange(main_figure[['ahc_t']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['ahc_t']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\ahc_t_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

pss_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\pSS\\pss_b.rds')
pss_b<-NormalizeData(pss_b)
pss_b<-FindVariableFeatures(pss_b,selection.method = 'vst',nfeatures = 2000)
pss_b<-ScaleData(pss_b)
pss_b<-RunPCA(pss_b,features = VariableFeatures(object = pss_b))
ElbowPlot(pss_b)
pss_b<-FindNeighbors(pss_b,dims = 1:15)
pss_b<-FindClusters(pss_b,resolution = 0.3)
pss_b<-RunUMAP(pss_b,dims = 1:15)
main_figure[['f6_j1']]<-DimPlot(pss_b,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(3,'Set1')))+
  labs(title = 'pSS PBMC B cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_j1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_J1.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

pss_b_expr_matrix<-as(as.matrix(pss_b@assays$RNA@counts),'sparseMatrix')
pss_b_p_data<-pss_b@meta.data
pss_b_f_data<-data.frame(gene_short_name=row.names(pss_b),row.names = row.names(pss_b))
pss_b_pd<-new('AnnotatedDataFrame',data=pss_b_p_data)
pss_b_fd<-new('AnnotatedDataFrame',data=pss_b_f_data)
pss_b_cds<-newCellDataSet(pss_b_expr_matrix,
                          phenoData = pss_b_pd,
                          featureData = pss_b_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
pss_b_cds<-estimateSizeFactors(pss_b_cds)
pss_b_cds<-estimateDispersions(pss_b_cds)
pss_b_cds<-detectGenes(pss_b_cds,min_expr = 0.5)
pss_b_disp_table<-dispersionTable(pss_b_cds)
pss_b_ordering_genes<-as.character(subset(pss_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
pss_b_cds<-setOrderingFilter(pss_b_cds,pss_b_ordering_genes)
pss_b_cds<-reduceDimension(pss_b_cds,max_components = 2,method = 'DDRTree')
pss_b_cds<-orderCells(pss_b_cds,reverse = F)
main_figure[['f6_k1']]<-plot_cell_trajectory(pss_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('pSS PBMC B cell subtypes trajectory')
main_figure[['f6_k1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_k2']]<-plot_cell_trajectory(pss_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('pSS PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f6_k2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

pss_b_time_diff<-differentialGeneTest(pss_b_cds[pss_b_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
pss_b_time_diff<-pss_b_time_diff[,c(5,2,3,4,1,6,7)]
pss_b_markers<-FindAllMarkers(pss_b,only.pos = TRUE,logfc.threshold = 0.5)
pss_b_top10<-pss_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
pss_b_top10_ordergene<-pss_b_time_diff[pss_b_top10$gene,]
pss_b_time_genes<-pss_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
pss_b_time_genes<-unique(pss_b_time_genes)
pss_b_time_genes<-na.omit(pss_b_time_genes)
pss_b_newdata<-data.frame(Pseudotime = seq(min(pData(pss_b_cds)$Pseudotime),
                                           max(pData(pss_b_cds)$Pseudotime),length.out=100))
pss_b_m<-genSmoothCurves(pss_b_cds[pss_b_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = pss_b_newdata)
pss_b_m<-pss_b_m[!apply(pss_b_m,1,sum)==0,]
pss_b_m<-log10(pss_b_m+1)
pss_b_m<-pss_b_m[!apply(pss_b_m,1,sd)==0,]
pss_b_m<-Matrix::t(scale(Matrix::t(pss_b_m),center = TRUE))
pss_b_m<-pss_b_m[is.na(row.names(pss_b_m))==FALSE,]
pss_b_m[is.nan(pss_b_m)]=0
pss_b_m[pss_b_m>3]=3
pss_b_m[pss_b_m<-3]=-3
pss_b_row_dist<-as.dist((1-cor(Matrix::t(pss_b_m)))/2)
pss_b_row_dist[is.na(pss_b_row_dist)]<-1

p1<-pheatmap(pss_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = pss_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

pss_b_annotation_col<-data.frame(pseudotime=rescale(pss_b_newdata$Pseudotime,to=c(-1,1)))
row.names(pss_b_annotation_col)<-colnames(pss_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
rownames(annotation_row)<-rownames(pss_b_m)
pss_b_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_l1']]<-pheatmap(pss_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = pss_b_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = pss_b_annotation_col,
                                 annotation_colors = pss_b_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_l1']]<-as.ggplot(main_figure[['f6_l1']])
main_figure[['f6_l1']]<-ggarrange(main_figure[['f6_l1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_l1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_L1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

sle_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\SLE\\sle_b.rds')
sle_b<-NormalizeData(sle_b)
sle_b<-FindVariableFeatures(sle_b,selection.method = 'vst',nfeatures = 2000)
sle_b<-ScaleData(sle_b)
sle_b<-RunPCA(sle_b,features = VariableFeatures(object = sle_b))
ElbowPlot(sle_b)
sle_b<-FindNeighbors(sle_b,dims = 1:15)
sle_b<-FindClusters(sle_b,resolution = 0.3)
sle_b<-RunUMAP(sle_b,dims = 1:15)
main_figure[['f6_j2']]<-DimPlot(sle_b,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(3,'Set1')))+
  labs(title = 'SLE PBMC B cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_j2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_J2.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

sle_b_expr_matrix<-as(as.matrix(sle_b@assays$RNA@counts),'sparseMatrix')
sle_b_p_data<-sle_b@meta.data
sle_b_f_data<-data.frame(gene_short_name=row.names(sle_b),row.names = row.names(sle_b))
sle_b_pd<-new('AnnotatedDataFrame',data=sle_b_p_data)
sle_b_fd<-new('AnnotatedDataFrame',data=sle_b_f_data)
sle_b_cds<-newCellDataSet(sle_b_expr_matrix,
                          phenoData = sle_b_pd,
                          featureData = sle_b_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
sle_b_cds<-estimateSizeFactors(sle_b_cds)
sle_b_cds<-estimateDispersions(sle_b_cds)
sle_b_cds<-detectGenes(sle_b_cds,min_expr = 0.5)
sle_b_disp_table<-dispersionTable(sle_b_cds)
sle_b_ordering_genes<-as.character(subset(sle_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
sle_b_cds<-setOrderingFilter(sle_b_cds,sle_b_ordering_genes)
sle_b_cds<-reduceDimension(sle_b_cds,max_components = 2,method = 'DDRTree')
sle_b_cds<-orderCells(sle_b_cds,reverse = F)
main_figure[['f6_k3']]<-plot_cell_trajectory(sle_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('SLE PBMC B cell subtypes trajectory')
main_figure[['f6_k3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K3.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_k4']]<-plot_cell_trajectory(sle_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('SLE PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f6_k4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K4.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

sle_b_time_diff<-differentialGeneTest(sle_b_cds[sle_b_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
sle_b_time_diff<-sle_b_time_diff[,c(5,2,3,4,1,6,7)]
sle_b_markers<-FindAllMarkers(sle_b,only.pos = TRUE,logfc.threshold = 0.5)
sle_b_top10<-sle_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
sle_b_top10_ordergene<-sle_b_time_diff[sle_b_top10$gene,]
sle_b_time_genes<-sle_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
sle_b_time_genes<-unique(sle_b_time_genes)
sle_b_time_genes<-na.omit(sle_b_time_genes)
sle_b_newdata<-data.frame(Pseudotime = seq(min(pData(sle_b_cds)$Pseudotime),
                                           max(pData(sle_b_cds)$Pseudotime),length.out=100))
sle_b_m<-genSmoothCurves(sle_b_cds[sle_b_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = sle_b_newdata)
sle_b_m<-sle_b_m[!apply(sle_b_m,1,sum)==0,]
sle_b_m<-log10(sle_b_m+1)
sle_b_m<-sle_b_m[!apply(sle_b_m,1,sd)==0,]
sle_b_m<-Matrix::t(scale(Matrix::t(sle_b_m),center = TRUE))
sle_b_m<-sle_b_m[is.na(row.names(sle_b_m))==FALSE,]
sle_b_m[is.nan(sle_b_m)]=0
sle_b_m[sle_b_m>3]=3
sle_b_m[sle_b_m<-3]=-3
sle_b_row_dist<-as.dist((1-cor(Matrix::t(sle_b_m)))/2)
sle_b_row_dist[is.na(sle_b_row_dist)]<-1

p1<-pheatmap(sle_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = sle_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

sle_b_annotation_col<-data.frame(pseudotime=rescale(sle_b_newdata$Pseudotime,to=c(-1,1)))
row.names(sle_b_annotation_col)<-colnames(sle_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(sle_b_m)
sle_b_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_l2']]<-pheatmap(sle_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = sle_b_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = sle_b_annotation_col,
                                 annotation_colors = sle_b_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_l2']]<-as.ggplot(main_figure[['f6_l2']])
main_figure[['f6_l2']]<-ggarrange(main_figure[['f6_l2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_l2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_L2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

ahc_b<-readRDS('D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_b.rds')
ahc_b<-NormalizeData(ahc_b)
ahc_b<-FindVariableFeatures(ahc_b,selection.method = 'vst',nfeatures = 2000)
ahc_b<-ScaleData(ahc_b)
ahc_b<-RunPCA(ahc_b,features = VariableFeatures(object = ahc_b))
ElbowPlot(ahc_b)
ahc_b<-FindNeighbors(ahc_b,dims = 1:15)
ahc_b<-FindClusters(ahc_b,resolution = 0.3)
ahc_b<-RunUMAP(ahc_b,dims = 1:15)
main_figure[['f6_j3']]<-DimPlot(ahc_b,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(3,'Set1')))+
  labs(title = 'aHC PBMC B cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_j3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_J3.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

ahc_b_expr_matrix<-as(as.matrix(ahc_b@assays$RNA@counts),'sparseMatrix')
ahc_b_p_data<-ahc_b@meta.data
ahc_b_f_data<-data.frame(gene_short_name=row.names(ahc_b),row.names = row.names(ahc_b))
ahc_b_pd<-new('AnnotatedDataFrame',data=ahc_b_p_data)
ahc_b_fd<-new('AnnotatedDataFrame',data=ahc_b_f_data)
ahc_b_cds<-newCellDataSet(ahc_b_expr_matrix,
                          phenoData = ahc_b_pd,
                          featureData = ahc_b_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
ahc_b_cds<-estimateSizeFactors(ahc_b_cds)
ahc_b_cds<-estimateDispersions(ahc_b_cds)
ahc_b_cds<-detectGenes(ahc_b_cds,min_expr = 0.5)
ahc_b_disp_table<-dispersionTable(ahc_b_cds)
ahc_b_ordering_genes<-as.character(subset(ahc_b_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
ahc_b_cds<-setOrderingFilter(ahc_b_cds,ahc_b_ordering_genes)
ahc_b_cds<-reduceDimension(ahc_b_cds,max_components = 2,method = 'DDRTree')
ahc_b_cds<-orderCells(ahc_b_cds,reverse = F)
main_figure[['f6_k5']]<-plot_cell_trajectory(ahc_b_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('aHC PBMC B cell subtypes trajectory')
main_figure[['f6_k5']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K5.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_k6']]<-plot_cell_trajectory(ahc_b_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('aHC PBMC B cell subtypes trajectory (Pseudotime)')
main_figure[['f6_k6']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_K6.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

ahc_b_time_diff<-differentialGeneTest(ahc_b_cds[ahc_b_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
ahc_b_time_diff<-ahc_b_time_diff[,c(5,2,3,4,1,6,7)]
ahc_b_markers<-FindAllMarkers(ahc_b,only.pos = TRUE,logfc.threshold = 0.5)
ahc_b_top10<-ahc_b_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
ahc_b_top10_ordergene<-ahc_b_time_diff[ahc_b_top10$gene,]
ahc_b_time_genes<-ahc_b_top10_ordergene %>% pull(gene_short_name) %>% as.character()
ahc_b_time_genes<-unique(ahc_b_time_genes)
ahc_b_time_genes<-na.omit(ahc_b_time_genes)
ahc_b_newdata<-data.frame(Pseudotime = seq(min(pData(ahc_b_cds)$Pseudotime),
                                           max(pData(ahc_b_cds)$Pseudotime),length.out=100))
ahc_b_m<-genSmoothCurves(ahc_b_cds[ahc_b_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = ahc_b_newdata)
ahc_b_m<-ahc_b_m[!apply(ahc_b_m,1,sum)==0,]
ahc_b_m<-log10(ahc_b_m+1)
ahc_b_m<-ahc_b_m[!apply(ahc_b_m,1,sd)==0,]
ahc_b_m<-Matrix::t(scale(Matrix::t(ahc_b_m),center = TRUE))
ahc_b_m<-ahc_b_m[is.na(row.names(ahc_b_m))==FALSE,]
ahc_b_m[is.nan(ahc_b_m)]=0
ahc_b_m[ahc_b_m>3]=3
ahc_b_m[ahc_b_m<-3]=-3
ahc_b_row_dist<-as.dist((1-cor(Matrix::t(ahc_b_m)))/2)
ahc_b_row_dist[is.na(ahc_b_row_dist)]<-1

p1<-pheatmap(ahc_b_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = ahc_b_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

ahc_b_annotation_col<-data.frame(pseudotime=rescale(ahc_b_newdata$Pseudotime,to=c(-1,1)))
row.names(ahc_b_annotation_col)<-colnames(ahc_b_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(ahc_b_m)
ahc_b_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['ahc_b']]<-pheatmap(ahc_b_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,clustering_method = 'ward.D2',
                                 clustering_distance_rows = ahc_b_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = ahc_b_annotation_col,
                                 annotation_colors = ahc_b_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['ahc_b']]<-as.ggplot(main_figure[['ahc_b']])
main_figure[['ahc_b']]<-ggarrange(main_figure[['ahc_b']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['ahc_b']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\ahc_b_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

pss_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\pSS\\pss_m.rds')
#pss_m<-subset(pss_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
pss_m<-NormalizeData(pss_m)
pss_m<-FindVariableFeatures(pss_m,selection.method = 'vst',nfeatures = 2000)
pss_m<-ScaleData(pss_m)
pss_m<-RunPCA(pss_m,features = VariableFeatures(object = pss_m))
ElbowPlot(pss_m)
pss_m<-FindNeighbors(pss_m,dims = 1:15)
pss_m<-FindClusters(pss_m,resolution = 0.3)
pss_m<-RunUMAP(pss_m,dims = 1:15)
main_figure[['f6_m1']]<-DimPlot(pss_m,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(4,'Set1')))+
  labs(title = 'pSS PBMC Myeloid cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_m1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_M1.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

pss_m_expr_matrix<-as(as.matrix(pss_m@assays$RNA@counts),'sparseMatrix')
pss_m_p_data<-pss_m@meta.data
pss_m_f_data<-data.frame(gene_short_name=row.names(pss_m),row.names = row.names(pss_m))
pss_m_pd<-new('AnnotatedDataFrame',data=pss_m_p_data)
pss_m_fd<-new('AnnotatedDataFrame',data=pss_m_f_data)
pss_m_cds<-newCellDataSet(pss_m_expr_matrix,
                          phenoData = pss_m_pd,
                          featureData = pss_m_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
pss_m_cds<-estimateSizeFactors(pss_m_cds)
pss_m_cds<-estimateDispersions(pss_m_cds)
pss_m_cds<-detectGenes(pss_m_cds,min_expr = 0.5)
pss_m_disp_table<-dispersionTable(pss_m_cds)
pss_m_ordering_genes<-as.character(subset(pss_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
pss_m_cds<-setOrderingFilter(pss_m_cds,pss_m_ordering_genes)
pss_m_cds<-reduceDimension(pss_m_cds,max_components = 2,method = 'DDRTree')
pss_m_cds<-orderCells(pss_m_cds,reverse = F)
main_figure[['f6_n1']]<-plot_cell_trajectory(pss_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('pSS PBMC Monocyte subtypes trajectory')
main_figure[['f6_n1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N1.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_n2']]<-plot_cell_trajectory(pss_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('pSS PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f6_n2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N2.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

pss_m_time_diff<-differentialGeneTest(pss_m_cds[pss_m_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
pss_m_time_diff<-pss_m_time_diff[,c(5,2,3,4,1,6,7)]
pss_m_markers<-FindAllMarkers(pss_m,only.pos = TRUE,logfc.threshold = 0.5)
pss_m_top10<-pss_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
pss_m_top10_ordergene<-pss_m_time_diff[pss_m_top10$gene,]
pss_m_time_genes<-pss_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
pss_m_time_genes<-unique(pss_m_time_genes)
pss_m_time_genes<-na.omit(pss_m_time_genes)
pss_m_newdata<-data.frame(Pseudotime = seq(min(pData(pss_m_cds)$Pseudotime),
                                           max(pData(pss_m_cds)$Pseudotime),length.out=100))
pss_m_m<-genSmoothCurves(pss_m_cds[pss_m_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = pss_m_newdata)
pss_m_m<-pss_m_m[!apply(pss_m_m,1,sum)==0,]
pss_m_m<-log10(pss_m_m+1)
pss_m_m<-pss_m_m[!apply(pss_m_m,1,sd)==0,]
pss_m_m<-Matrix::t(scale(Matrix::t(pss_m_m),center = TRUE))
pss_m_m<-pss_m_m[is.na(row.names(pss_m_m))==FALSE,]
pss_m_m[is.nan(pss_m_m)]=0
pss_m_m[pss_m_m>3]=3
pss_m_m[pss_m_m<-3]=-3
pss_m_row_dist<-as.dist((1-cor(Matrix::t(pss_m_m)))/2)
pss_m_row_dist[is.na(pss_m_row_dist)]<-1

p1<-pheatmap(pss_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = pss_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

pss_m_annotation_col<-data.frame(pseudotime=rescale(pss_m_newdata$Pseudotime,to=c(-1,1)))
row.names(pss_m_annotation_col)<-colnames(pss_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(pss_m_m)
pss_m_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_o1']]<-pheatmap(pss_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = pss_m_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = pss_m_annotation_col,
                                 annotation_colors = pss_m_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_o1']]<-as.ggplot(main_figure[['f6_o1']])
main_figure[['f6_o1']]<-ggarrange(main_figure[['f6_o1']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_o1']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_O1.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

sle_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\SLE\\sle_m.rds')
#sle_m<-subset(sle_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
sle_m<-NormalizeData(sle_m)
sle_m<-FindVariableFeatures(sle_m,selection.method = 'vst',nfeatures = 2000)
sle_m<-ScaleData(sle_m)
sle_m<-RunPCA(sle_m,features = VariableFeatures(object = sle_m))
ElbowPlot(sle_m)
sle_m<-FindNeighbors(sle_m,dims = 1:15)
sle_m<-FindClusters(sle_m,resolution = 0.3)
sle_m<-RunUMAP(sle_m,dims = 1:15)
main_figure[['f6_m2']]<-DimPlot(sle_m,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(4,'Set1')))+
  labs(title = 'SLE PBMC Myeloid cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_m2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_M2.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

sle_m_expr_matrix<-as(as.matrix(sle_m@assays$RNA@counts),'sparseMatrix')
sle_m_p_data<-sle_m@meta.data
sle_m_f_data<-data.frame(gene_short_name=row.names(sle_m),row.names = row.names(sle_m))
sle_m_pd<-new('AnnotatedDataFrame',data=sle_m_p_data)
sle_m_fd<-new('AnnotatedDataFrame',data=sle_m_f_data)
sle_m_cds<-newCellDataSet(sle_m_expr_matrix,
                          phenoData = sle_m_pd,
                          featureData = sle_m_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
sle_m_cds<-estimateSizeFactors(sle_m_cds)
sle_m_cds<-estimateDispersions(sle_m_cds)
sle_m_cds<-detectGenes(sle_m_cds,min_expr = 0.5)
sle_m_disp_table<-dispersionTable(sle_m_cds)
sle_m_ordering_genes<-as.character(subset(sle_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
sle_m_cds<-setOrderingFilter(sle_m_cds,sle_m_ordering_genes)
sle_m_cds<-reduceDimension(sle_m_cds,max_components = 2,method = 'DDRTree')
sle_m_cds<-orderCells(sle_m_cds,reverse = F)
main_figure[['f6_n3']]<-plot_cell_trajectory(sle_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('SLE PBMC Monocyte subtypes trajectory')
main_figure[['f6_n3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N3.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_n4']]<-plot_cell_trajectory(sle_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('SLE PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f6_n4']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N4.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")


sle_m_time_diff<-differentialGeneTest(sle_m_cds[sle_m_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
sle_m_time_diff<-sle_m_time_diff[,c(5,2,3,4,1,6,7)]
sle_m_markers<-FindAllMarkers(sle_m,only.pos = TRUE,logfc.threshold = 0.5)
sle_m_top10<-sle_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
sle_m_top10_ordergene<-sle_m_time_diff[sle_m_top10$gene,]
sle_m_time_genes<-sle_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
sle_m_time_genes<-unique(sle_m_time_genes)
sle_m_time_genes<-na.omit(sle_m_time_genes)
sle_m_newdata<-data.frame(Pseudotime = seq(min(pData(sle_m_cds)$Pseudotime),
                                           max(pData(sle_m_cds)$Pseudotime),length.out=100))
sle_m_m<-genSmoothCurves(sle_m_cds[sle_m_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = sle_m_newdata)
sle_m_m<-sle_m_m[!apply(sle_m_m,1,sum)==0,]
sle_m_m<-log10(sle_m_m+1)
sle_m_m<-sle_m_m[!apply(sle_m_m,1,sd)==0,]
sle_m_m<-Matrix::t(scale(Matrix::t(sle_m_m),center = TRUE))
sle_m_m<-sle_m_m[is.na(row.names(sle_m_m))==FALSE,]
sle_m_m[is.nan(sle_m_m)]=0
sle_m_m[sle_m_m>3]=3
sle_m_m[sle_m_m<-3]=-3
sle_m_row_dist<-as.dist((1-cor(Matrix::t(sle_m_m)))/2)
sle_m_row_dist[is.na(sle_m_row_dist)]<-1

p1<-pheatmap(sle_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = sle_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

sle_m_annotation_col<-data.frame(pseudotime=rescale(sle_m_newdata$Pseudotime,to=c(-1,1)))
row.names(sle_m_annotation_col)<-colnames(sle_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(sle_m_m)
sle_m_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['f6_o2']]<-pheatmap(sle_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,clustering_method = 'ward.D2',
                                 show_colnames = FALSE,cutree_rows = 4,
                                 clustering_distance_rows = sle_m_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = sle_m_annotation_col,
                                 annotation_colors = sle_m_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['f6_o2']]<-as.ggplot(main_figure[['f6_o2']])
main_figure[['f6_o2']]<-ggarrange(main_figure[['f6_o2']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['f6_o2']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_O2.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")

ahc_m<-readRDS('D:\\scRNA\\JIA\\trajectory\\aHC\\ahc_m.rds')
ahc_m<-NormalizeData(ahc_m)
ahc_m<-FindVariableFeatures(ahc_m,selection.method = 'vst',nfeatures = 2000)
ahc_m<-ScaleData(ahc_m)
ahc_m<-RunPCA(ahc_m,features = VariableFeatures(object = ahc_m))
ElbowPlot(ahc_m)
ahc_m<-FindNeighbors(ahc_m,dims = 1:15)
ahc_m<-FindClusters(ahc_m,resolution = 0.3)
ahc_m<-RunUMAP(ahc_m,dims = 1:15)
main_figure[['f6_m3']]<-DimPlot(ahc_m,group.by = 'celltype_2',label = TRUE,cols = c(brewer.pal(4,'Set1')))+
  labs(title = 'aHC PBMC Myeloid cell subtypes',x='UMAP1',y='UMAP2')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA,color = 'black'),
        aspect.ratio = 1)
main_figure[['f6_m3']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_M3.pdf',units = 'in',dpi = 300,
       width = 6,height = 6,device = "pdf")

ahc_m<-subset(ahc_m,subset = celltype_2 %in% c('CD14 Mono','Inter Mono','CD16 Mono'))
ahc_m_expr_matrix<-as(as.matrix(ahc_m@assays$RNA@counts),'sparseMatrix')
ahc_m_p_data<-ahc_m@meta.data
ahc_m_f_data<-data.frame(gene_short_name=row.names(ahc_m),row.names = row.names(ahc_m))
ahc_m_pd<-new('AnnotatedDataFrame',data=ahc_m_p_data)
ahc_m_fd<-new('AnnotatedDataFrame',data=ahc_m_f_data)
ahc_m_cds<-newCellDataSet(ahc_m_expr_matrix,
                          phenoData = ahc_m_pd,
                          featureData = ahc_m_fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
ahc_m_cds<-estimateSizeFactors(ahc_m_cds)
ahc_m_cds<-estimateDispersions(ahc_m_cds)
ahc_m_cds<-detectGenes(ahc_m_cds,min_expr = 0.5)
ahc_m_disp_table<-dispersionTable(ahc_m_cds)
ahc_m_ordering_genes<-as.character(subset(ahc_m_disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
ahc_m_cds<-setOrderingFilter(ahc_m_cds,ahc_m_ordering_genes)
ahc_m_cds<-reduceDimension(ahc_m_cds,max_components = 2,method = 'DDRTree')
ahc_m_cds<-orderCells(ahc_m_cds,reverse = F)
main_figure[['f6_n5']]<-plot_cell_trajectory(ahc_m_cds,show_cell_names = FALSE,color_by = 'celltype_2')+
  ggtitle('aHC PBMC Monocyte subtypes trajectory')
main_figure[['f6_n5']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N5.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

main_figure[['f6_n6']]<-plot_cell_trajectory(ahc_m_cds,show_cell_names = FALSE,color_by = 'Pseudotime')+
  ggtitle('aHC PBMC Monocyte subtypes trajectory (Pseudotime)')
main_figure[['f6_n6']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\F6_N6.pdf',units = 'in',dpi = 300,
       width = 9,height = 6,device = "pdf")

ahc_m_time_diff<-differentialGeneTest(ahc_m_cds[ahc_m_ordering_genes,],cores = 1,
                                      fullModelFormulaStr = '~sm.ns(Pseudotime)')
ahc_m_time_diff<-ahc_m_time_diff[,c(5,2,3,4,1,6,7)]
ahc_m_markers<-FindAllMarkers(ahc_m,only.pos = TRUE,logfc.threshold = 0.5)
ahc_m_top10<-ahc_m_markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
ahc_m_top10_ordergene<-ahc_m_time_diff[ahc_m_top10$gene,]
ahc_m_time_genes<-ahc_m_top10_ordergene %>% pull(gene_short_name) %>% as.character()
ahc_m_time_genes<-unique(ahc_m_time_genes)
ahc_m_time_genes<-na.omit(ahc_m_time_genes)
ahc_m_newdata<-data.frame(Pseudotime = seq(min(pData(ahc_m_cds)$Pseudotime),
                                           max(pData(ahc_m_cds)$Pseudotime),length.out=100))
ahc_m_m<-genSmoothCurves(ahc_m_cds[ahc_m_time_genes],
                         trend_formula = '~sm.ns(Pseudotime,df=3)',
                         relative_expr = T,new_data = ahc_m_newdata)
ahc_m_m<-ahc_m_m[!apply(ahc_m_m,1,sum)==0,]
ahc_m_m<-log10(ahc_m_m+1)
ahc_m_m<-ahc_m_m[!apply(ahc_m_m,1,sd)==0,]
ahc_m_m<-Matrix::t(scale(Matrix::t(ahc_m_m),center = TRUE))
ahc_m_m<-ahc_m_m[is.na(row.names(ahc_m_m))==FALSE,]
ahc_m_m[is.nan(ahc_m_m)]=0
ahc_m_m[ahc_m_m>3]=3
ahc_m_m[ahc_m_m<-3]=-3
ahc_m_row_dist<-as.dist((1-cor(Matrix::t(ahc_m_m)))/2)
ahc_m_row_dist[is.na(ahc_m_row_dist)]<-1

p1<-pheatmap(ahc_m_m,useRater = TRUE,cluster_cols = FALSE,cluster_rows = TRUE,
             show_rownames = FALSE,show_colnames = FALSE,clustering_method = 'ward.D2',
             clustering_distance_rows = ahc_m_row_dist,cutree_rows = 4,
             border_color = NA,filename = NA,
             color = colorRampPalette(c("navy","white","firebrick3"))(100))

ahc_m_annotation_col<-data.frame(pseudotime=rescale(ahc_m_newdata$Pseudotime,to=c(-1,1)))
row.names(ahc_m_annotation_col)<-colnames(ahc_m_m)
annotation_row<-data.frame(Cluster=factor(cutree(p1$tree_row,4)))
row.names(annotation_row)<-rownames(ahc_m_m)
ahc_m_anno_colors<-list(pseudotime=viridis(100),
                        Cluster=row_color)
main_figure[['ahc_m']]<-pheatmap(ahc_m_m,useRaster = T,cluster_rows = TRUE,
                                 cluster_cols = FALSE,show_rownames = TRUE,
                                 show_colnames = FALSE,cutree_rows = 4,clustering_method = 'ward.D2',
                                 clustering_distance_rows = ahc_m_row_dist,
                                 border_color = NA,filename = NA,
                                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                                 annotation_col = ahc_m_annotation_col,
                                 annotation_colors = ahc_m_anno_colors,
                                 annotation_row = annotation_row)
main_figure[['ahc_m']]<-as.ggplot(main_figure[['ahc_m']])
main_figure[['ahc_m']]<-ggarrange(main_figure[['ahc_m']],ncol = 1,nrow = 1,common.legend = TRUE,legend = 'right')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = 'cm'))
main_figure[['ahc_m']]
ggsave('D:\\scRNA\\JIA\\JIA主figure-6-11\\ahc_m_heatmap.pdf',units = 'in',dpi = 300,
       width = 9,height = 12,device = "pdf")
