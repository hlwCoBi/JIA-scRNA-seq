de_doublet<-function(seurat_object,batch_ident,PCuse=1:30,exp_doublet_ratio=0.07){
  Idents(seurat_object)<-batch_ident
  result<-data.frame(cellname=NA,doublet_result=NA)
  result<-result[0,]
  for(i in levels(Idents(seurat_object))){
    sub_seu<-subset(seurat_object,idents = i)
    sub_seu<-NormalizeData(sub_seu)
    sub_seu<-FindVariableFeatures(sub_seu)
    sub_seu<-ScaleData(sub_seu)
    sub_seu<-RunPCA(sub_seu,dims=PCuse)
    sub_seu<-RunUMAP(sub_seu,dims=PCuse)
    sub_seu<-FindNeighbors(sub_seu,dims=PCuse)
    sub_seu<-FindClusters(sub_seu,resolution = 0.4)
    sweep.res.test<-paramSweep_v3(sub_seu,PCs = PCuse,sct=F)
    sweep.stats.test<-summarizeSweep(sweep.res.test,GT=F)
    bcmvn.test<-find.pK(sweep.stats.test)
    bcmvn.test<-dplyr::arrange(bcmvn.test,desc(BCmetric))
    best_pK<-as.numeric(bcmvn.test[1,"pK"])
    homotypic.prop<-modelHomotypic(sub_seu@active.ident)
    nExp_poi <- round(exp_doublet_ratio*ncol(sub_seu))
    nExp_poi.adj<-round(nExp_poi*(1-homotypic.prop))
    sub_seu<-doubletFinder_v3(seu = sub_seu,PCs = PCuse,pN = 0.25,pK = best_pK,nExp = nExp_poi.adj)
    sub_result<-data.frame(cellname=colnames(sub_seu),doublet_result=sub_seu@meta.data[,ncol(sub_seu@meta.data)])
    result<-rbind(result,sub_result)
  }
  rownames(result)<-result$cellname
  seurat_object$doublet<-result[colnames(seurat_object),"doublet_result"]
  seurat_object@meta.data$doublet<-factor(seurat_object$doublet,levels = c('Singlet','Doublet'))
  return(seurat_object)
}
