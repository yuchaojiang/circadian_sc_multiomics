# count: a vector of expression or accessibility across cells
# cell.meta: the name of the metainfo of cells from the sc object
aggregateCell=function(count, sc, cell.meta, FUN, FUN.name=NULL, ...){
  if(is.null(FUN.name)){FUN.name=FUN}
  if(!all(cell.meta %in% colnames(sc@meta.data))){
    stop('All cell.meta elements need to be in sc@meta.data!')
  }
  if(length(count)!=ncol(sc)) {
    stop('Length of count needs to be the same as the number 
         of cells in the single-cell object!')
  }
  output=aggregate(count, by=sc@meta.data[,cell.meta], FUN=FUN)
  colnames(output)[ncol(output)]=FUN.name
  # output=output[complete.cases(output),] # Remove any NA, NAN etc. Not here to make the genes aligned.
  return(output)
}


runJTK=function(output.list, group.by, ...){
  if(!group.by %in% colnames(output.list[[1]])){
    stop('group.by must be a column name of output.list[[1]]!')
  }
  output=matrix(nrow=length(output.list), ncol=nrow(output.list[[1]]))
  rownames(output)=names(output.list)
  colnames(output)=apply(output.list[[1]][,-ncol(output.list[[1]])],1, function(x){paste(x, collapse = '_')})
  for(i in 1:nrow(output)){
    output[i,]=output.list[[i]][,ncol(output.list[[i]])]
  }
  output.meta=output.list[[1]][,-ncol(output.list[[1]])]
  
  JTK=vector(mode = "list", length=length(unique(output.meta[,group.by])))
  names(JTK)=unique(output.meta[,group.by])
  
  for(i in 1:length(JTK)){
    group.i=names(JTK)[i]
    output.i=output[,output.meta[,group.by]==group.i]
    output.meta.i=output.meta[output.meta[,group.by]==group.i,]
    output.i=data.frame(gene=names(output.list),output.i)
    
    write.csv(output.i, file='circadian/temp.csv',row.names = FALSE)
    JTK.i=meta2d(infile = 'circadian/temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = output.meta.i$ZT, cycMethod = c("JTK"),
                 maxper = 24, minper = 8, combinePvalue = 'bonferroni')$meta
    file.remove('circadian/temp.csv')
    JTK[[i]]=cbind(group=group.i, JTK.i)
  }
  return(JTK)
}



runJTK.FUN=function(output.list, group.by, cell.meta, FUN, ...){
  if(!group.by %in% colnames(output.list[[1]])){
    stop('group.by must be a column name of output.list[[1]]!')
  }
  output=matrix(nrow=length(output.list), ncol=nrow(output.list[[1]]))
  rownames(output)=names(output.list)
  colnames(output)=apply(output.list[[1]][,cell.meta],1, function(x){paste(x, collapse = '_')})
  for(i in 1:nrow(output)){
    output[i,]=output.list[[i]][,FUN]
  }
  output.meta=output.list[[1]][,cell.meta]
  
  JTK=vector(mode = "list", length=length(unique(output.meta[,group.by])))
  names(JTK)=unique(output.meta[,group.by])
  
  for(i in 1:length(JTK)){
    group.i=names(JTK)[i]
    output.i=output[,output.meta[,group.by]==group.i]
    if(FUN=='cv' | FUN=='nonzero.mean'){
      output.i[is.nan(output.i)]=0
    }
    output.meta.i=output.meta[output.meta[,group.by]==group.i,]
    output.i=data.frame(gene=names(output.list),output.i)
    
    write.csv(output.i, file='circadian/temp.csv',row.names = FALSE)
    JTK.i=meta2d(infile = 'circadian/temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = output.meta.i$ZT, cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni')$meta
    file.remove('circadian/temp.csv')
    JTK[[i]]=cbind(group=group.i, JTK.i)
  }
  return(JTK)
}



keep_major_celltype=function(sc, from=NULL, to=NULL){
  if(is.null(from)){
    from='seurat_clusters'
  }
  if(is.null(to)){
    to='celltype'
  }
  temp=table(as.matrix(sc@meta.data[from]), as.matrix(sc@meta.data[to]))
  keep=c()
  for(i in 1:nrow(temp)){
    keep=c(keep, paste(rownames(temp)[i], colnames(temp)[which.max(temp[i,])]))
  }
  sc=sc[,paste(as.matrix(sc@meta.data[from]), as.matrix(sc@meta.data[to])) %in% keep]
  return(sc)
}

rna_count_for_jtk_input = function(sc){
  tic()
  sc@assays$RNA@counts %>% as.data.frame() %>% 
    mutate(across(everything(), 
                  function(x){(x/sum(x))*median(sc$nCount_RNA)})) -> df_
  
  genes_ = c()
  sc@meta.data[,c("celltype", "ZT", "cluster")] %>% 
    as.data.frame() %>% rownames_to_column(var = "cell") %>% 
    group_by(celltype, ZT, cluster) %>% 
    group_map(function(x,y){
      cbind(y,x) -> df_meta
      celltype = unique(df_meta$celltype)
      ZT = unique(df_meta$ZT)
      cluster = unique(df_meta$cluster)
      #    df_[, df_meta$cell][1:5,1:5] %>% 
      df_[df_meta$cell] %>% 
        apply(., 1, mean) %>% 
        as.data.frame() %>% "colnames<-"(., "mean") %>% 
        mutate(celltype = celltype, ZT = ZT, cluster = cluster) %>% 
        rownames_to_column("gene") %>% 
        dplyr::select(gene, celltype, ZT, cluster, mean) 
    }) %>% do.call(rbind, .) %>% group_by(gene) %>% 
    group_map(function(x,y){
      c(genes_, y$gene) ->> genes_
      x %>% as.data.frame()
    }) -> output.rna.list
  names(output.rna.list) <- genes_
  rm(df_, genes_)
  toc()
  return(output.rna.list)
}

### Don't directly run this function
### Wraping furrr::map into a function cause ram outrage.
### copy the code to R script and run
atac_count_for_jtk_input = function(sc, peak.ref, gene.promoter.ref){
  tic()
  findOverlaps(peak.ref, gene.promoter.ref) %>% 
    as.data.frame() -> df_overlaps
  
  gene.promoter.ref %>% data.frame() %>% .[, "symbol"] %>% 
    as.data.frame() %>% 
    "colnames<-"(., "symbol") %>% 
    rownames_to_column("subjectIndex") %>% 
    mutate(subjectIndex = as.integer(subjectIndex)) -> df_index
  left_join(x = df_overlaps, y = df_index, by = c("subjectHits" = "subjectIndex")) -> df_index
  
  sc@assays$ATAC@counts -> atac_count
  sc$nCount_ATAC -> sc_nCount_ATAC
  sc@meta.data -> sc_meta
#  rownames(sc@assays$RNA@counts) -> genes
  
  3500*1024^2 -> max_size
  options(future.globals.maxSize= max_size)
  plan(multisession, workers = 6)
  
  furrr::future_map(unique(df_index$subjectHits), function(x){
    df_index[df_index$subjectHits == x, ] -> x
    x$queryHits -> idxs
    unique(x$symbol) -> symbol
    
    atac_count[idxs, , drop = FALSE] %>% 
      #  sc@assays$ATAC@counts[idxs, , drop = FALSE] %>% 
      #      colSums() %>% {(./sc$nCount_ATAC[1:10])*median(sc$nCount_ATAC)} -> exp
      colSums() %>% {(./sc_nCount_ATAC)*median(sc_nCount_ATAC)} -> exp
    sc_meta[names(exp), c("celltype", "ZT", "cluster")] -> df_meta
    df_meta %>% mutate(exp = exp) %>% 
      group_by(celltype, ZT, cluster) %>% 
      summarise(mean = mean(exp)) %>% as.data.frame() -> df_
    attr(x = df_, which = "gene") <- symbol
    return(df_)
  }) -> output.atac.list
  map(output.atac.list, function(x){
    attributes(x)$gene
  }) %>% unlist() -> names(output.atac.list)
  
#  setdiff(genes, names(output.atac.list)) -> rest_genes
#  output.atac.list[[1]] %>% mutate(mean = 0) -> empty_DF
#  map(1:length(rest_genes), function(x){
#    empty_DF
#  }) -> rest_genes_list
#  names(rest_genes_list) <- rest_genes
#  c(output.atac.list, rest_genes_list) -> output.atac.list
  
  toc()
  return(output.atac.list)
}

gene_activity_for_jtk_input = function(sc){
  tic()
  sc@assays$gene_activity@counts %>% as.data.frame() %>% 
    mutate(across(everything(), 
                  function(x){(x/sum(x))*median(sc$nCount_gene_activity)})) -> df_
  
  genes_ = c()
  sc@meta.data[,c("celltype", "ZT", "cluster")] %>% 
    as.data.frame() %>% rownames_to_column(var = "cell") %>% 
    group_by(celltype, ZT, cluster) %>% 
    group_map(function(x,y){
      cbind(y,x) -> df_meta
      celltype = unique(df_meta$celltype)
      ZT = unique(df_meta$ZT)
      cluster = unique(df_meta$cluster)
      #    df_[, df_meta$cell][1:5,1:5] %>% 
      df_[, df_meta$cell] %>% 
        apply(., 1, mean) %>% 
        as.data.frame() %>% "colnames<-"(., "mean") %>% 
        mutate(celltype = celltype, ZT = ZT, cluster = cluster) %>% 
        rownames_to_column("gene") %>% 
        dplyr::select(gene, celltype, ZT, cluster, mean) 
    }) %>% do.call(rbind, .) %>% group_by(gene) %>% 
    group_map(function(x,y){
      c(genes_, y$gene) ->> genes_
      x %>% as.data.frame()
    }) -> output.gene.activity.list
  names(output.gene.activity.list) <- genes_
  rm(df_, genes_)
  toc()
  return(output.gene.activity.list)
}

TF_score_for_jtk_input = function(sc){
  tic()
  sc@assays$MOTIF@data %>% as.data.frame() -> df_
  
  genes_ = c()
  sc@meta.data[,c("celltype", "ZT", "cluster")] %>% 
    as.data.frame() %>% rownames_to_column(var = "cell") %>% 
    group_by(celltype, ZT, cluster) %>% 
    group_map(function(x,y){
      cbind(y,x) -> df_meta
      celltype = unique(df_meta$celltype)
      ZT = unique(df_meta$ZT)
      cluster = unique(df_meta$cluster)
      #    df_[, df_meta$cell][1:5,1:5] %>% 
      df_[, df_meta$cell] %>% 
        apply(., 1, mean) %>% 
        as.data.frame() %>% "colnames<-"(., "mean") %>% 
        mutate(celltype = celltype, ZT = ZT, cluster = cluster) %>% 
        rownames_to_column("gene") %>% 
        dplyr::select(gene, celltype, ZT, cluster, mean) 
    }) %>% do.call(rbind, .) %>% group_by(gene) %>% 
    group_map(function(x,y){
      c(genes_, y$gene) ->> genes_
      x %>% as.data.frame()
    }) -> output.TF.score.list
  names(output.TF.score.list) <- genes_
  rm(df_, genes_)
  toc()
  return(output.TF.score.list)
}

gene_activity_for_jtk_input_1 = function(sc, assay_name){
  tic()
  print("start generating normalized df")
  sc@assays[[assay_name]]@counts %>% as.data.frame() %>% 
    mutate(across(everything(), 
                  function(x){(x/sum(x))*median(sc[[sprintf("nCount_%s", assay_name)]] %>% unlist())})) -> df_
  print("normalized df generated")
  
  genes_ = c()
  sc@meta.data[,c("celltype", "ZT", "cluster")] %>% 
    as.data.frame() %>% rownames_to_column(var = "cell") %>% 
    group_by(celltype, ZT, cluster) %>% 
    group_map(function(x,y){
      cbind(y,x) -> df_meta
      celltype = unique(df_meta$celltype)
      ZT = unique(df_meta$ZT)
      cluster = unique(df_meta$cluster)
      print(sprintf("Processing: celltype=%s, ZT=%s, cluster=%s", 
                    celltype, ZT, cluster))
      #    df_[, df_meta$cell][1:5,1:5] %>% 
      df_[, df_meta$cell] %>% 
        apply(., 1, mean) %>% 
        as.data.frame() %>% "colnames<-"(., "mean") %>% 
        mutate(celltype = celltype, ZT = ZT, cluster = cluster) %>% 
        rownames_to_column("gene") %>% 
        dplyr::select(gene, celltype, ZT, cluster, mean) 
    }) %>% do.call(rbind, .) %>% group_by(gene) %>% 
    group_map(function(x,y){
      c(genes_, y$gene) ->> genes_
      x %>% as.data.frame()
    }) -> output.gene.activity.list
  names(output.gene.activity.list) <- genes_
  rm(df_, genes_)
  toc()
  return(output.gene.activity.list)
}

moment_fun<-function(Y.temp,cellsize=NULL){
  if(is.null(cellsize)){cellsize<-rep(1, length(Y.temp))}
  m1<-sum(Y.temp)/sum(cellsize)
  m2<-sum(Y.temp*(Y.temp-1))/sum(cellsize^2)
  m3<-sum(Y.temp*(Y.temp-1)*(Y.temp-2))/sum(cellsize^3)
  kon.hat<--2*(-m1*m2^2+m1^2*m3)/(-m1*m2^2+2*m1^2*m3-m2*m3)
  koff.hat<-2*(m1^2-m2)*(m1*m2-m3)*(m1*m3-m2^2)/(m1^2*m2-2*m2^2+m1*m3)/(2*m1^2*m3-m1*m2^2-m2*m3)
  s.hat<-(-m1*m2^2+2*m1^2*m3-m2*m3)/(m1^2*m2-2*m2^2+m1*m3)
  kinetic.estimates<-list(kon.hat=round(kon.hat,4),
                          koff.hat=round(koff.hat,4),
                          scale=round(s.hat,2))
  return(kinetic.estimates)
}


moment_fun_kon<-function(Y.temp,cellsize=NULL){
  if(is.null(cellsize)){cellsize<-rep(1, length(Y.temp))}
  m1<-sum(Y.temp)/sum(cellsize)
  m2<-sum(Y.temp*(Y.temp-1))/sum(cellsize^2)
  m3<-sum(Y.temp*(Y.temp-1)*(Y.temp-2))/sum(cellsize^3)
  kon.hat<--2*(-m1*m2^2+m1^2*m3)/(-m1*m2^2+2*m1^2*m3-m2*m3)
  if(is.nan(kon.hat)){kon.hat=0}
  return(kon.hat)
}


