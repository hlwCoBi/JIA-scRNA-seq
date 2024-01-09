tcells_jia<-subset(tcells,subset = group %in% c('b27 negative','b27 positive','cHC'))
tcells_jia<-NormalizeData(tcells_jia)
tcells_jia<-FindVariableFeatures(tcells_jia,selection.method = 'vst',nfeatures = 2000)
tcells_jia<-ScaleData(tcells_jia)
tcells_jia<-RunPCA(tcells_jia,features = VariableFeatures(object = tcells_jia))
ElbowPlot(tcells_jia)
tcells_jia<-FindNeighbors(tcells_jia,dims = 1:20)
tcells_jia<-FindClusters(tcells_jia,resolution = 0.3)
tcells_jia<-RunUMAP(tcells_jia,dims = 1:20)
DimPlot_plus(tcells_jia,group.by = NULL)
DimPlot(tcells_jia,group.by = 'group',label = TRUE)
DotPlot(tcells_jia,features = c('NKG7','KLRF1','CD3E','CD4','CD40LG','CD8A','CCR7','SELL',
                                'PTPRC','IL7R','TCF7','LEF1','FOXP3','IL2RA','TRAV1-2',
                                'CD27','CD28','CD44','KLRG1','TRDC'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:14,celltype='NA')
tcells_jia_markers<-FindAllMarkers(tcells_jia,only.pos = TRUE)

for (i in 0:14) {
  t<-subset(tcells_jia_markers,cluster == i)
  write.csv(t,paste0('D:\\scRNA\\JIA\\JIA-9-20\\tcells\\tcells_jia\\markers\\cluster_',i,'.csv'))
}
celltype_2<-data.frame(clusterid=0:14,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,4),2]<-'CD4 TCM'
celltype_2[celltype_2$clusterid %in% c(1,7),2]<-'CD4 Naive'
celltype_2[celltype_2$clusterid %in% c(2,13),2]<-'CCR7+ T'
celltype_2[celltype_2$clusterid %in% c(3),2]<-'NK'
celltype_2[celltype_2$clusterid %in% c(5,9,12),2]<-'CD8 Naive'
celltype_2[celltype_2$clusterid %in% c(6,8,10,11),2]<-'CD8 TEM'
for (i in 1:nrow(celltype_2)) {
  tcells_jia@meta.data[which(tcells_jia@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}

tcells_SLE<-subset(tcells,subset = group == 'SLE')
tcells_SLE<-NormalizeData(tcells_SLE)
tcells_SLE<-FindVariableFeatures(tcells_SLE,selection.method = 'vst',nfeatures = 2000)
tcells_SLE<-ScaleData(tcells_SLE)
tcells_SLE<-RunPCA(tcells_SLE,features = VariableFeatures(object = tcells_SLE))
ElbowPlot(tcells_SLE)
tcells_SLE<-FindNeighbors(tcells_SLE,dims = 1:20)
tcells_SLE<-FindClusters(tcells_SLE,resolution = 0.3)
tcells_SLE<-RunUMAP(tcells_SLE,dims = 1:20)
DimPlot_plus(tcells_SLE,group.by = NULL)
DotPlot(tcells_SLE,features = c('NKG7','KLRF1','CD3E','CD4','CD40LG','CD8A','CCR7','SELL',
                                'PTPRC','IL7R','TCF7','LEF1','FOXP3','IL2RA','TRAV1-2',
                                'CD27','CD28','CD44','KLRG1','TRDC'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:7,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,6),2]<-'CD8 TEM'
celltype_2[celltype_2$clusterid %in% c(1,4),2]<-'CD8 Naive'
celltype_2[celltype_2$clusterid %in% c(2),2]<-'CD4 Naive'
celltype_2[celltype_2$clusterid %in% c(3,5),2]<-'NK'
celltype_2[celltype_2$clusterid %in% c(7),2]<-'CD4 TCM'
for (i in 1:nrow(celltype_2)) {
  tcells_SLE@meta.data[which(tcells_SLE@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(tcells_SLE@meta.data$celltype_2)) {
  c<-subset(tcells_SLE,subset = celltype_2 == i)
  for (j in 1:nrow(tcells@meta.data)) {
    if(rownames(tcells@meta.data)[j] %in% rownames(c@meta.data)){
      tcells@meta.data$celltype_2[j]<-i
    }
  }
}

tcells_pSS<-subset(tcells,subset = group == 'pSS')
tcells_pSS<-NormalizeData(tcells_pSS)
tcells_pSS<-FindVariableFeatures(tcells_pSS,selection.method = 'vst',nfeatures = 2000)
tcells_pSS<-ScaleData(tcells_pSS)
tcells_pSS<-RunPCA(tcells_pSS,features = VariableFeatures(object = tcells_pSS))
ElbowPlot(tcells_pSS)
tcells_pSS<-FindNeighbors(tcells_pSS,dims = 1:20)
tcells_pSS<-FindClusters(tcells_pSS,resolution = 0.3)
tcells_pSS<-RunUMAP(tcells_pSS,dims = 1:20)
DimPlot_plus(tcells_pSS,group.by = NULL)
DotPlot(tcells_pSS,features = c('NKG7','KLRF1','CD3E','CD4','CD40LG','CD8A','CCR7','SELL',
                                'PTPRC','IL7R','TCF7','LEF1','FOXP3','IL2RA','TRAV1-2',
                                'CD27','CD28','CD44','KLRG1','TRDC'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:15,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,3,12,14,15),2]<-'NK'
celltype_2[celltype_2$clusterid %in% c(1),2]<-'CD4 TCM'
celltype_2[celltype_2$clusterid %in% c(2,13),2]<-'CD8 TEM'
celltype_2[celltype_2$clusterid %in% c(4,6,11),2]<-'CD4 Naive'
celltype_2[celltype_2$clusterid %in% c(5),2]<-'CD8 Naive'
celltype_2[celltype_2$clusterid %in% c(7),2]<-'CD4 TEM'
celltype_2[celltype_2$clusterid %in% c(8),2]<-'MAIT'
celltype_2[celltype_2$clusterid %in% c(9),2]<-'Treg'
celltype_2[celltype_2$clusterid %in% c(10),2]<-'ILC'

for (i in 1:nrow(celltype_2)) {
  tcells_pSS@meta.data[which(tcells_pSS@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(tcells_pSS@meta.data$celltype_2)) {
  c<-subset(tcells_pSS,subset = celltype_2 == i)
  for (j in 1:nrow(tcells@meta.data)) {
    if(rownames(tcells@meta.data)[j] %in% rownames(c@meta.data)){
      tcells@meta.data$celltype_2[j]<-i
    }
  }
}

tcells_aHC<-subset(tcells,subset = group == 'aHC')
tcells_aHC<-NormalizeData(tcells_aHC)
tcells_aHC<-FindVariableFeatures(tcells_aHC,selection.method = 'vst',nfeatures = 2000)
tcells_aHC<-ScaleData(tcells_aHC)
tcells_aHC<-RunPCA(tcells_aHC,features = VariableFeatures(object = tcells_aHC))
ElbowPlot(tcells_aHC)
tcells_aHC<-FindNeighbors(tcells_aHC,dims = 1:20)
tcells_aHC<-FindClusters(tcells_aHC,resolution = 0.3)
tcells_aHC<-RunUMAP(tcells_aHC,dims = 1:20)
DimPlot_plus(tcells_aHC,group.by = NULL)
DotPlot(tcells_aHC,features = c('NKG7','KLRF1','CD3E','CD4','CD40LG','CD8A','CCR7','SELL',
                                'PTPRC','IL7R','TCF7','LEF1','FOXP3','IL2RA','TRAV1-2',
                                'CD27','CD28','CD44','KLRG1','TRDC'),group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:11,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0),2]<-'CD4 Naive'
celltype_2[celltype_2$clusterid %in% c(1),2]<-'CD4 TCM'
celltype_2[celltype_2$clusterid %in% c(2,8),2]<-'CD8 Naive'
celltype_2[celltype_2$clusterid %in% c(3,11),2]<-'CD8 TEM'
celltype_2[celltype_2$clusterid %in% c(4,5,10),2]<-'NK'
celltype_2[celltype_2$clusterid %in% c(6),2]<-'ILC'
celltype_2[celltype_2$clusterid %in% c(7),2]<-'MAIT'
celltype_2[celltype_2$clusterid %in% c(9),2]<-'PTPRC+ T'
for (i in 1:nrow(celltype_2)) {
  tcells_aHC@meta.data[which(tcells_aHC@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(tcells_aHC@meta.data$celltype_2)) {
  c<-subset(tcells_aHC,subset = celltype_2 == i)
  for (j in 1:nrow(tcells@meta.data)) {
    if(rownames(tcells@meta.data)[j] %in% rownames(c@meta.data)){
      tcells@meta.data$celltype_2[j]<-i
    }
  }
}

bcells<-subset(seurat,subset = celltype_1 == 'B cells')
bcells<-NormalizeData(bcells)
bcells<-FindVariableFeatures(bcells,selection.method = 'vst',nfeatures = 2000)
bcells<-ScaleData(bcells)
bcells<-RunPCA(bcells,features = VariableFeatures(object = bcells))
ElbowPlot(bcells)
bcells<-FindNeighbors(bcells,dims = 1:20)
bcells<-FindClusters(bcells,resolution = 0.3)
bcells<-RunUMAP(bcells,dims = 1:20)
DimPlot_plus(bcells,group.by = NULL)
DimPlot(bcells,group.by = 'group',label = TRUE)
bcells_markers<-FindAllMarkers(bcells,only.pos = TRUE)
for (i in 0:13) {
  b<-subset(bcells_markers,cluster==i)
  write.csv(b,paste0('D:\\scRNA\\JIA\\JIA-9-20\\bcells\\markers\\cluster_',i,'.csv'))
}
DotPlot(bcells,features = c('CD38','MZB1','XBP1','CD79A','CD79B','CD19','MS4A1','IGHD',
                            'TCL1A','CD27','TNFRSF13B'),
        cols = c('white','red'),group.by = NULL)+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))

bcells_jia<-subset(bcells,subset = group %in% c('b27 negative','b27 positive','cHC'))
bcells_jia<-NormalizeData(bcells_jia)
bcells_jia<-FindVariableFeatures(bcells_jia,selection.method = 'vst',nfeatures = 2000)
bcells_jia<-ScaleData(bcells_jia)
bcells_jia<-RunPCA(bcells_jia,features = VariableFeatures(object = bcells_jia))
ElbowPlot(bcells_jia)
bcells_jia<-FindNeighbors(bcells_jia,dims = 1:20)
bcells_jia<-FindClusters(bcells_jia,resolution = 0.3)
bcells_jia<-RunUMAP(bcells_jia,dims = 1:20)
bcells_jia_markers<-FindAllMarkers(bcells_jia,only.pos = TRUE)
DimPlot_plus(bcells_jia,group.by = NULL)
DimPlot(bcells_jia,group.by = 'group',label = TRUE)
DotPlot(bcells_jia,features = c('CD38','MZB1','XBP1','CD27','TNFRSF13B','CD79A','CD79B',
                                'CD19','MS4A1','IGHD','TCL1A'),
        cols = c('white','red'),group.by = NULL)+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:11,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,1,4,5,6,8,10,11),2]<-'Naive B'
celltype_2[celltype_2$clusterid %in% c(2,3),2]<-'Memory B'
celltype_2[celltype_2$clusterid %in% c(7,9),2]<-'Plasma B'
for (i in 1:nrow(celltype_2)) {
  bcells_jia@meta.data[which(bcells_jia@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
bcells@meta.data$celltype_2<-'NA'

for (i in unique(bcells_jia@meta.data$celltype_2)) {
  c<-subset(bcells_jia,subset = celltype_2 == i)
  for (j in 1:nrow(bcells@meta.data)) {
    if(rownames(bcells@meta.data)[j] %in% rownames(c@meta.data)){
      bcells@meta.data$celltype_2[j]<-i
    }
  }
}

bcells_sle<-subset(bcells,subset = group == 'SLE')
bcells_sle<-NormalizeData(bcells_sle)
bcells_sle<-FindVariableFeatures(bcells_sle,selection.method = 'vst',nfeatures = 2000)
bcells_sle<-ScaleData(bcells_sle)
bcells_sle<-RunPCA(bcells_sle,features = VariableFeatures(object = bcells_sle))
ElbowPlot(bcells_sle)
bcells_sle<-FindNeighbors(bcells_sle,dims = 1:20)
bcells_sle<-FindClusters(bcells_sle,resolution = 0.3)
bcells_sle<-RunUMAP(bcells_sle,dims = 1:20)
bcells_sle_markers<-FindAllMarkers(bcells_sle,only.pos = TRUE)
for (i in 0:7) {
  b<-subset(bcells_sle_markers,cluster==i)
  write.csv(b,paste0('D:\\scRNA\\JIA\\JIA-9-20\\bcells\\bcells_sle\\markers\\cluster_',i,'.csv'))
}
DimPlot_plus(bcells_sle,group.by = NULL)
DotPlot(bcells_sle,features = c('CD38','MZB1','XBP1','CD27','TNFRSF13B','CD79A','CD79B',
                                'CD19','MS4A1','IGHD','TCL1A'),
        cols = c('white','red'),group.by = NULL)+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:7,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,2,3,5),2]<-'Naive B'
celltype_2[celltype_2$clusterid %in% c(1,6),2]<-'Memory B'
celltype_2[celltype_2$clusterid %in% c(4,7),2]<-'Plasma B'
for (i in 1:nrow(celltype_2)) {
  bcells_sle@meta.data[which(bcells_sle@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(bcells_sle@meta.data$celltype_2)) {
  c<-subset(bcells_sle,subset = celltype_2 == i)
  for (j in 1:nrow(bcells@meta.data)) {
    if(rownames(bcells@meta.data)[j] %in% rownames(c@meta.data)){
      bcells@meta.data$celltype_2[j]<-i
    }
  }
}

bcells_pss<-subset(bcells,subset = group == 'pSS')
bcells_pss<-NormalizeData(bcells_pss)
bcells_pss<-FindVariableFeatures(bcells_pss,selection.method = 'vst',nfeatures = 2000)
bcells_pss<-ScaleData(bcells_pss)
bcells_pss<-RunPCA(bcells_pss,features = VariableFeatures(object = bcells_pss))
ElbowPlot(bcells_pss)
bcells_pss<-FindNeighbors(bcells_pss,dims = 1:20)
bcells_pss<-FindClusters(bcells_pss,resolution = 0.3)
bcells_pss<-RunUMAP(bcells_pss,dims = 1:20)
bcells_pss_markers<-FindAllMarkers(bcells_pss,only.pos = TRUE)
for (i in 0:5) {
  b<-subset(bcells_pss_markers,cluster==i)
  write.csv(b,paste0('D:\\scRNA\\JIA\\JIA-9-20\\bcells\\bcells_pss\\markers\\cluster_',i,'.csv'))
}
DimPlot_plus(bcells_pss,group.by = NULL)
DotPlot(bcells_pss,features = c('CD38','MZB1','XBP1','CD27','TNFRSF13B','CD79A','CD79B',
                                'CD19','MS4A1','IGHD','TCL1A'),
        cols = c('white','red'),group.by = NULL)+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:5,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,1,3),2]<-'Naive B'
celltype_2[celltype_2$clusterid %in% c(2),2]<-'Memory B'
celltype_2[celltype_2$clusterid %in% c(4,5),2]<-'Plasma B'
for (i in 1:nrow(celltype_2)) {
  bcells_pss@meta.data[which(bcells_pss@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(bcells_pss@meta.data$celltype_2)) {
  c<-subset(bcells_pss,subset = celltype_2 == i)
  for (j in 1:nrow(bcells@meta.data)) {
    if(rownames(bcells@meta.data)[j] %in% rownames(c@meta.data)){
      bcells@meta.data$celltype_2[j]<-i
    }
  }
}

bcells_ahc<-subset(bcells,subset = group == 'aHC')
bcells_ahc<-NormalizeData(bcells_ahc)
bcells_ahc<-FindVariableFeatures(bcells_ahc,selection.method = 'vst',nfeatures = 2000)
bcells_ahc<-ScaleData(bcells_ahc)
bcells_ahc<-RunPCA(bcells_ahc,features = VariableFeatures(object = bcells_ahc))
ElbowPlot(bcells_ahc)
bcells_ahc<-FindNeighbors(bcells_ahc,dims = 1:20)
bcells_ahc<-FindClusters(bcells_ahc,resolution = 0.3)
bcells_ahc<-RunUMAP(bcells_ahc,dims = 1:20)
bcells_ahc_markers<-FindAllMarkers(bcells_ahc,only.pos = TRUE)
for (i in 0:4) {
  b<-subset(bcells_ahc_markers,cluster==i)
  write.csv(b,paste0('D:\\scRNA\\JIA\\JIA-9-20\\bcells\\bcells_ahc\\markers\\cluster_',i,'.csv'))
}
DimPlot_plus(bcells_ahc,group.by = NULL)
DotPlot(bcells_ahc,features = c('CD38','MZB1','XBP1','CD27','TNFRSF13B','CD79A','CD79B',
                                'CD19','MS4A1','IGHD','TCL1A'),
        cols = c('white','red'),group.by = NULL)+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:4,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,3),2]<-'Naive B'
celltype_2[celltype_2$clusterid %in% c(1),2]<-'Memory B'
celltype_2[celltype_2$clusterid %in% c(2,4),2]<-'Plasma B'
for (i in 1:nrow(celltype_2)) {
  bcells_ahc@meta.data[which(bcells_ahc@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
for (i in unique(bcells_ahc@meta.data$celltype_2)) {
  c<-subset(bcells_ahc,subset = celltype_2 == i)
  for (j in 1:nrow(bcells@meta.data)) {
    if(rownames(bcells@meta.data)[j] %in% rownames(c@meta.data)){
      bcells@meta.data$celltype_2[j]<-i
    }
  }
}

myeloidcells<-subset(seurat,subset = celltype_1 == 'Myeloid cells')
myeloidcells<-NormalizeData(myeloidcells)
myeloidcells<-FindVariableFeatures(myeloidcells,selection.method = 'vst',nfeatures = 2000)
myeloidcells<-ScaleData(myeloidcells)
myeloidcells<-RunPCA(myeloidcells,features = VariableFeatures(object = myeloidcells))
ElbowPlot(myeloidcells)
myeloidcells<-FindNeighbors(myeloidcells,dims = 1:20)
myeloidcells<-FindClusters(myeloidcells,resolution = 0.3)
myeloidcells<-RunUMAP(myeloidcells,dims = 1:20)
DimPlot_plus(myeloidcells,group.by = NULL)
DimPlot(myeloidcells,group.by = 'group',label = TRUE)
myeloidcells_markers<-FindAllMarkers(myeloidcells,only.pos = TRUE)
for (i in 0:18) {
  m<-subset(myeloidcells_markers,cluster==i)
  write.csv(m,paste0('D:\\scRNA\\JIA\\JIA-9-20\\myeloidcells\\markers\\cluster_',i,'.csv'))
}
DotPlot(myeloidcells,features = c('CD68','MS4A7','CD4','LYZ','CD14','FCGR3A','FCER1A','CD1C',
                                  'CLEC4C'),group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))

myeloidcells_jia<-subset(myeloidcells,
                         subset = group %in% c('b27 negative','b27 positive','cHC'))
myeloidcells_jia<-NormalizeData(myeloidcells_jia)
myeloidcells_jia<-FindVariableFeatures(myeloidcells_jia,selection.method = 'vst',nfeatures = 2000)
myeloidcells_jia<-ScaleData(myeloidcells_jia)
myeloidcells_jia<-RunPCA(myeloidcells_jia,
                         features = VariableFeatures(object = myeloidcells_jia))
ElbowPlot(myeloidcells_jia)
myeloidcells_jia<-FindNeighbors(myeloidcells_jia,dims = 1:20)
myeloidcells_jia<-FindClusters(myeloidcells_jia,resolution = 0.3)
myeloidcells_jia<-RunUMAP(myeloidcells_jia,dims = 1:20)
myeloidcells_jia_markers<-FindAllMarkers(myeloidcells_jia,only.pos = TRUE)
for (i in 0:14) {
  m<-subset(myeloidcells_jia_markers,cluster==i)
  write.csv(m,paste0('D:\\scRNA\\JIA\\JIA-9-20\\myeloidcells\\myeloidcells_jia\\markers\\cluster',i,'.csv'))
}
DimPlot_plus(myeloidcells_jia,group.by = NULL)
DimPlot(myeloidcells_jia,group.by = 'group',label = TRUE)
DotPlot(myeloidcells_jia,features = c('CD68','MS4A7','CD4','LYZ','CD14','FCGR3A','FCER1A',
                                      'CST3','CD1C','CLEC4C','CXCL2','CX3CR1','IL1B','PLTP',
                                      'NLRP3'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))

myeloidcells@meta.data$celltype_2<-'NA'
celltype_2<-data.frame(clusterid=0:14,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,1,2,4,5,9,10,12),2]<-'CD14 Mono'
celltype_2[celltype_2$clusterid %in% c(3,13,14),2]<-'Inter Mono'
celltype_2[celltype_2$clusterid %in% c(6),2]<-'CD14 Mono'
celltype_2[celltype_2$clusterid %in% c(7),2]<-'DC'
celltype_2[celltype_2$clusterid %in% c(8,11),2]<-'CD16 Mono'
for (i in 1:nrow(celltype_2)) {
  myeloidcells_jia@meta.data[which(myeloidcells_jia@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}

myeloidcells_pss<-subset(myeloidcells,subset = group == 'pSS')
myeloidcells_pss<-NormalizeData(myeloidcells_pss)
myeloidcells_pss<-FindVariableFeatures(myeloidcells_pss,selection.method = 'vst',nfeatures = 2000)
myeloidcells_pss<-ScaleData(myeloidcells_pss)
myeloidcells_pss<-RunPCA(myeloidcells_pss,features = VariableFeatures(object = myeloidcells_pss))
ElbowPlot(myeloidcells_pss)
myeloidcells_pss<-FindNeighbors(myeloidcells_pss,dims = 1:20)
myeloidcells_pss<-FindClusters(myeloidcells_pss,resolution = 0.3)
myeloidcells_pss<-RunUMAP(myeloidcells_pss,dims = 1:20)
myeloidcells_pss_markers<-FindAllMarkers(myeloidcells_pss,only.pos = TRUE)
for (i in 0:8) {
  m<-subset(myeloidcells_pss_markers,cluster==i)
  write.csv(m,paste0('D:\\scRNA\\JIA\\JIA-9-20\\myeloidcells\\myeloidcells_pss\\markers\\cluster_',i,'.csv'))
}
DimPlot_plus(myeloidcells_pss,group.by = NULL)
DotPlot(myeloidcells_pss,features = c('CD68','MS4A7','CD4','LYZ','CD14','FCGR3A','FCER1A',
                                      'CST3','CD1C','CLEC4C','CXCL2','CX3CR1','IL1B','PLTP',
                                      'NLRP3'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:8,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,1,3,7,8),2]<-'CD14 Mono'
celltype_2[celltype_2$clusterid %in% c(2,5),2]<-'CD16 Mono'
celltype_2[celltype_2$clusterid %in% c(4),2]<-'Inter Mono'
celltype_2[celltype_2$clusterid %in% c(6),2]<-'DC'
for (i in 1:nrow(celltype_2)) {
  myeloidcells_pss@meta.data[which(myeloidcells_pss@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}

myeloidcells_sle<-subset(myeloidcells,subset = group == 'SLE')
myeloidcells_sle<-NormalizeData(myeloidcells_sle)
myeloidcells_sle<-FindVariableFeatures(myeloidcells_sle,selection.method = 'vst',nfeatures = 2000)
myeloidcells_sle<-ScaleData(myeloidcells_sle)
myeloidcells_sle<-RunPCA(myeloidcells_sle,features = VariableFeatures(object = myeloidcells_sle))
ElbowPlot(myeloidcells_sle)
myeloidcells_sle<-FindNeighbors(myeloidcells_sle,dims = 1:20)
myeloidcells_sle<-FindClusters(myeloidcells_sle,resolution = 0.3)
myeloidcells_sle<-RunUMAP(myeloidcells_sle,dims = 1:20)
myeloidcells_sle_markers<-FindAllMarkers(myeloidcells_sle,only.pos = TRUE)
for (i in 0:11) {
  m<-subset(myeloidcells_sle_markers,cluster==i)
  write.csv(m,paste0('D:\\scRNA\\JIA\\JIA-9-20\\myeloidcells\\myeloidcells_sle\\markers\\cluster_',i,'.csv'))
}
DimPlot_plus(myeloidcells_sle,group.by = NULL)
DotPlot(myeloidcells_sle,features = c('CD68','MS4A7','CD4','LYZ','CD14','FCGR3A','FCER1A',
                                      'CST3','CD1C','CLEC4C','CXCL2','CX3CR1','IL1B','PLTP',
                                      'NLRP3'),group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:11,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,2,3,4,5),2]<-'CD14 Mono'
celltype_2[celltype_2$clusterid %in% c(1,7,9),2]<-'CD16 Mono'
celltype_2[celltype_2$clusterid %in% c(6),2]<-'DC'
celltype_2[celltype_2$clusterid %in% c(8,10,11),2]<-'Inter Mono'
for (i in 1:nrow(celltype_2)) {
  myeloidcells_sle@meta.data[which(myeloidcells_sle@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}

myeloidcells_aHC<-subset(myeloidcells,subset = group == 'aHC')
myeloidcells_aHC<-NormalizeData(myeloidcells_aHC)
myeloidcells_aHC<-FindVariableFeatures(myeloidcells_aHC,selection.method = 'vst',nfeatures = 2000)
myeloidcells_aHC<-ScaleData(myeloidcells_aHC)
myeloidcells_aHC<-RunPCA(myeloidcells_aHC,features = VariableFeatures(object = myeloidcells_aHC))
ElbowPlot(myeloidcells_aHC)
myeloidcells_aHC<-FindNeighbors(myeloidcells_aHC,dims = 1:20)
myeloidcells_aHC<-FindClusters(myeloidcells_aHC,resolution = 0.3)
myeloidcells_aHC<-RunUMAP(myeloidcells_aHC,dims = 1:20)
DimPlot_plus(myeloidcells_aHC,group.by = NULL)
DotPlot(myeloidcells_aHC,features = c('CD68','MS4A7','CD4','LYZ','CD14','FCGR3A','FCER1A',
                                      'CST3','CD1C','CLEC4C','CXCL2','CX3CR1','IL1B','PLTP',
                                      'NLRP3'),
        group.by = NULL,cols = c('white','red'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
celltype_2<-data.frame(clusterid=0:4,celltype='NA')
celltype_2[celltype_2$clusterid %in% c(0,1),2]<-'CD14 Mono'
celltype_2[celltype_2$clusterid %in% c(2),2]<-'Inter Mono'
celltype_2[celltype_2$clusterid %in% c(3),2]<-'CD16 Mono'
celltype_2[celltype_2$clusterid %in% c(4),2]<-'DC'
for (i in 1:nrow(celltype_2)) {
  myeloidcells_aHC@meta.data[which(myeloidcells_aHC@meta.data$seurat_clusters==celltype_2$clusterid[i]),'celltype_2']<-celltype_2$celltype[i]
}
