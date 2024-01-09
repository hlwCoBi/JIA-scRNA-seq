matrix_gene_filter<-function(expr_matrix,quantile=0.03,min_expr=5,return_gene=T){
  min_sum<-ncol(expr_matrix)*quantile*min_expr
  pass_gene<-rownames(expr_matrix)[Matrix::rowSums(expr_matrix)>=min_sum]
  expr_matrix<-expr_matrix[Matrix::rowSums(expr_matrix)>=min_sum,]
  quantile_res<-apply(expr_matrix,MARGIN = 1,FUN = stats::quantile,probs = (1-quantile),na.rm = F)
  if(!return_gene){
    return(expr_matrix[quantile_res>min_expr,])
  }else{
      return(rownames(expr_matrix)[quantile_res>min_expr])
    }
}

expr_rename_10x<-function(expr_matrix,prefix){
  meta<-data.frame(old_name=colnames(expr_matrix),new_name=NA,stringsAsFactors = F)
  meta$new_name<-paste0(prefix,1:ncol(expr_matrix))
  colnames(expr_matrix)<-meta$new_name
  return(list(expr=expr_matrix,meta=meta))
}

umi_gene_count<-function(expr_matrix){
  umi_count<-Matrix::colSums(expr_matrix)
  gene_count<-(expr_matrix>0)*1
  gene_count<-Matrix::colSums(gene_count)
  return(data.frame(umi_count=umi_count,gene_count=gene_count))
}
  
