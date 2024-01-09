seurat<-CreateSeuratObject(counts = all_expr,object='JIA',min.cells = 0,min.features = 0)
seurat@meta.data$group<-all_meta$group
seurat@meta.data$person<-all_meta$person
seurat@meta.data$umi_count<-all_meta$umi_count
seurat@meta.data$gene_count<-all_meta$gene_count
seurat<-subset(seurat,subset = gene_count>1000 & umi_count>1500 & percent.mt<10)
seurat<-de_doublet(seurat_object = seurat,batch_ident = seurat@meta.data$person,PCuse = 1:30,exp_doublet_ratio = 0.07)
seurat<-subset(seurat,subset = doublet == 'Singlet')
seurat<-seurat[-grep(pattern = '^MT-',x=rownames(seurat)),]
seurat<-NormalizeData(seurat)
seurat<-FindVariableFeatures(seurat,selection.method = 'vst',nfeatures = 2000)
seurat<-ScaleData(seurat)
seurat<-RunPCA(seurat,features = VariableFeatures(object = seurat))
ElbowPlot(seurat)
seurat<-FindNeighbors(seurat,dims = 1:20)
seurat<-FindClusters(seurat,resolution = 0.3)
seurat<-RunUMAP(seurat,dims = 1:20)
all_markers<-FindAllMarkers(seurat,only.pos = TRUE)
for (i in 0:26) {
  c<-all_markers[all_markers$cluster==i,]
  write.csv(c,file = paste0("all_markers/cluster_",i,".csv"))
}
celltype_1<-data.frame(clusterid=0:26,celltype='NA')
celltype_1[celltype_1$clusterid %in% c(0,3,5,7,11),2]<-'T cells'
celltype_1[celltype_1$clusterid %in% c(6,8,9,10,12,22),2]<-'NK cells'
celltype_1[celltype_1$clusterid %in% c(2,18,20,21,24,25),2]<-'B cells'
celltype_1[celltype_1$clusterid %in% c(1,4,13,14,15,16,17,19,23),2]<-'Myeloid cells'
celltype_1[celltype_1$clusterid %in% c(26),2]<-'Neutrophil'
for (i in 1:nrow(celltype_1)) {
  seurat@meta.data[which(seurat@meta.data$seurat_clusters==celltype_1$clusterid[i]),'celltype_1']<-celltype_1$celltype[i]
}
test<-CellSelector(plot = DimPlot(subset(seurat,seurat_clusters==7)))
for (i in 1:nrow(seurat@meta.data)) {
  if(rownames(seurat@meta.data)[i] %in% test){
    seurat@meta.data$celltype_1[i]<-'Platelet'
  }
}
